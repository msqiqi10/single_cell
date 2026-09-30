#!/usr/bin/env bash
set -euo pipefail

ROOT="/home/zzz0054/bio3"
PY="$ROOT/.venv/bin/python"
BUILD="$ROOT/notebooks/scripts/iNKT/build_inkt_legacy_ppt_qc_notebook.py"
EXECUTE="$ROOT/notebooks/scripts/iNKT/execute_inkt_notebook.py"
VALIDATE_QC="$ROOT/notebooks/scripts/iNKT/validate_inkt_legacy_ppt_qc.py"
SIG="$ROOT/notebooks/scripts/iNKT/analyze_inkt_gene_signatures.py"
SIG_TEST="$ROOT/notebooks/scripts/iNKT/test_analyze_inkt_gene_signatures.py"
DE="$ROOT/notebooks/scripts/iNKT/run_inkt_de_pathway_overlap.py"
TRAJ="$ROOT/notebooks/scripts/iNKT/analyze_inkt_trajectory_biomarkers.py"
CLUSTER_NAMESPACE="legacy_ppt_qc_leiden_res_0_5"

if [[ "${1:-}" != "--inside-tmux" ]]; then
    if [[ ! -x "$PY" ]]; then
        echo "Missing project Python: $PY" >&2
        exit 2
    fi
    if ! command -v tmux >/dev/null 2>&1; then
        echo "tmux is required for this long-running pipeline" >&2
        exit 2
    fi

    stamp="$(date -u +%Y%m%d_%H%M%S)"
    session="inkt_legacy_qc_${stamp}"
    run_root="$ROOT/output/iNKT_legacy_ppt_qc_runs/$stamp"
    preprocess_dir="$run_root/preprocess"
    notebook="$run_root/notebooks/scanpy_iNKT_legacy_ppt_qc_pipeline.ipynb"
    log="$run_root/tmux/pipeline.log"

    mkdir -p "$run_root/tmux" "$preprocess_dir"
    "$PY" "$BUILD" \
        --output-notebook "$notebook" \
        --run-dir "$preprocess_dir"

    tmux new-session -d \
        -s "$session" \
        -n pipeline \
        "env INKT_LEGACY_RUN_ROOT='$run_root' INKT_LEGACY_NOTEBOOK='$notebook' INKT_LEGACY_SESSION='$session' bash '$ROOT/scripts/runners/run_inkt_legacy_ppt_qc_pipeline.sh' --inside-tmux >'$log' 2>&1"

    printf '%s\n' "$session" > "$ROOT/.inkt_legacy_ppt_qc_session"
    printf '%s\n' "$run_root" > "$ROOT/.inkt_legacy_ppt_qc_run_dir"

    echo "Started legacy-PPT-QC pipeline."
    echo "Session: $session"
    echo "Run root: $run_root"
    echo "Log: $log"
    exit 0
fi

: "${INKT_LEGACY_RUN_ROOT:?INKT_LEGACY_RUN_ROOT is required}"
: "${INKT_LEGACY_NOTEBOOK:?INKT_LEGACY_NOTEBOOK is required}"
: "${INKT_LEGACY_SESSION:?INKT_LEGACY_SESSION is required}"

RUN_ROOT="$INKT_LEGACY_RUN_ROOT"
NOTEBOOK="$INKT_LEGACY_NOTEBOOK"
PREPROCESS_DIR="$RUN_ROOT/preprocess"
H5AD="$PREPROCESS_DIR/inkt_scanpy_tutorial_processed.h5ad"
EXTENDED_DIR="$RUN_ROOT/extended"
LOG_DIR="$RUN_ROOT/tmux/logs"
STATUS_DIR="$RUN_ROOT/tmux/status"

mkdir -p "$EXTENDED_DIR" "$LOG_DIR" "$STATUS_DIR"
cd "$ROOT"

export CUDA_VISIBLE_DEVICES=""
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
export MKL_NUM_THREADS=4
export NUMEXPR_NUM_THREADS=4
export NUMBA_NUM_THREADS=4

pipeline_started="$(date -u +%FT%TZ)"
printf '%s\n' "$pipeline_started" > "$RUN_ROOT/started_utc.txt"

finish_pipeline() {
    rc=$?
    printf '%s\n' "$rc" > "$STATUS_DIR/pipeline.exit"
    date -u +%FT%TZ > "$RUN_ROOT/finished_utc.txt"
    exit "$rc"
}
trap finish_pipeline EXIT

echo "[1/5] Executing preprocessing, clustering, plots, and notebook trajectory"
"$PY" "$EXECUTE" --notebook "$NOTEBOOK" --timeout 7200

echo "[2/5] Validating exact legacy-PPT QC reconstruction"
"$PY" "$VALIDATE_QC" \
    --h5ad "$H5AD" \
    --filter-summary "$PREPROCESS_DIR/tables/filter_summary.json" \
    --output "$RUN_ROOT/qc_validation.json"

echo "[3/5] Running code and input preflight checks"
"$PY" -m py_compile \
    "$BUILD" \
    "$EXECUTE" \
    "$VALIDATE_QC" \
    "$SIG" \
    "$SIG_TEST" \
    "$DE" \
    "$TRAJ"
"$PY" "$SIG_TEST" -v > "$LOG_DIR/signature_tests.log" 2>&1
"$PY" "$SIG" --h5ad "$H5AD" --validate-only > "$LOG_DIR/signature_validate.log" 2>&1
"$PY" "$TRAJ" \
    --input-h5ad "$H5AD" \
    --output-dir "$EXTENDED_DIR/trajectory" \
    --cluster-namespace "$CLUSTER_NAMESPACE" \
    --validate-only > "$LOG_DIR/trajectory_validate.log" 2>&1
"$PY" "$DE" --help > "$LOG_DIR/de_help.log" 2>&1

declare -A pids
jobs=(gene_sets signature_scores standardized_legacy marker_validation legacy de0 de1 trajectory)

start_job() {
    name="$1"
    threads="$2"
    shift 2
    (
        set +e
        env \
            CUDA_VISIBLE_DEVICES="" \
            OMP_NUM_THREADS="$threads" \
            OPENBLAS_NUM_THREADS="$threads" \
            MKL_NUM_THREADS="$threads" \
            NUMEXPR_NUM_THREADS="$threads" \
            NUMBA_NUM_THREADS="$threads" \
            "$@" > "$LOG_DIR/${name}.log" 2>&1
        rc=$?
        printf '%s\n' "$rc" > "$STATUS_DIR/${name}.exit"
        exit "$rc"
    ) &
    pids["$name"]=$!
}

echo "[4/5] Running signatures, marker validation, DE shards, legacy extraction, and trajectory"
start_job gene_sets 1 \
    "$PY" "$SIG" \
    --task gene-sets --threads 1 \
    --h5ad "$H5AD" \
    --out-dir "$EXTENDED_DIR/signatures"

start_job signature_scores 3 \
    "$PY" "$SIG" \
    --task signature-scores --threads 3 \
    --h5ad "$H5AD" \
    --out-dir "$EXTENDED_DIR/signatures"

start_job standardized_legacy 3 \
    "$PY" "$SIG" \
    --task standardized-legacy --threads 3 --bootstrap 10000 \
    --h5ad "$H5AD" \
    --out-dir "$EXTENDED_DIR/signatures"

start_job marker_validation 2 \
    "$PY" "$SIG" \
    --task marker-validation --threads 2 \
    --h5ad "$H5AD" \
    --out-dir "$EXTENDED_DIR/signatures"

start_job legacy 1 \
    "$PY" "$DE" \
    --task legacy \
    --h5ad "$H5AD" \
    --cluster-namespace "$CLUSTER_NAMESPACE" \
    --outdir "$EXTENDED_DIR/de_pathway"

start_job de0 1 \
    "$PY" "$DE" \
    --task current-de --n-jobs 7 \
    --shard-count 2 --shard-index 0 \
    --h5ad "$H5AD" \
    --cluster-namespace "$CLUSTER_NAMESPACE" \
    --outdir "$EXTENDED_DIR/de_pathway"

start_job de1 1 \
    "$PY" "$DE" \
    --task current-de --n-jobs 7 \
    --shard-count 2 --shard-index 1 \
    --h5ad "$H5AD" \
    --cluster-namespace "$CLUSTER_NAMESPACE" \
    --outdir "$EXTENDED_DIR/de_pathway"

start_job trajectory 3 \
    "$PY" "$TRAJ" \
    --input-h5ad "$H5AD" \
    --cluster-namespace "$CLUSTER_NAMESPACE" \
    --output-dir "$EXTENDED_DIR/trajectory"

overall=0
for name in "${jobs[@]}"; do
    if ! wait "${pids[$name]}"; then
        overall=1
    fi
done

dependency_ok=1
for name in legacy de0 de1; do
    if [[ ! -f "$STATUS_DIR/${name}.exit" ]] || [[ "$(<"$STATUS_DIR/${name}.exit")" != "0" ]]; then
        dependency_ok=0
    fi
done

echo "[5/5] Integrating overlap and offline enrichment"
if [[ "$dependency_ok" -eq 1 ]]; then
    set +e
    env \
        CUDA_VISIBLE_DEVICES="" \
        OMP_NUM_THREADS=1 \
        OPENBLAS_NUM_THREADS=1 \
        MKL_NUM_THREADS=1 \
        NUMEXPR_NUM_THREADS=1 \
        NUMBA_NUM_THREADS=1 \
        "$PY" "$DE" \
        --task overlap-pathway \
        --h5ad "$H5AD" \
        --cluster-namespace "$CLUSTER_NAMESPACE" \
        --outdir "$EXTENDED_DIR/de_pathway" \
        > "$LOG_DIR/overlap_pathway.log" 2>&1
    overlap_rc=$?
    set -e
else
    overlap_rc=90
    echo "Skipped because legacy or DE shard failed" > "$LOG_DIR/overlap_pathway.log"
fi
printf '%s\n' "$overlap_rc" > "$STATUS_DIR/overlap_pathway.exit"
if [[ "$overlap_rc" -ne 0 ]]; then
    overall=1
fi

{
    printf 'started_utc=%s\n' "$pipeline_started"
    date -u +completed_utc=%FT%TZ
    printf 'qc_validation=0\n'
    for name in "${jobs[@]}" overlap_pathway; do
        if [[ -f "$STATUS_DIR/${name}.exit" ]]; then
            printf '%s=%s\n' "$name" "$(<"$STATUS_DIR/${name}.exit")"
        else
            printf '%s=missing\n' "$name"
            overall=1
        fi
    done
    printf 'overall=%s\n' "$overall"
} > "$RUN_ROOT/pipeline_status.txt"

if [[ "$overall" -ne 0 ]]; then
    echo "Pipeline failed; inspect $RUN_ROOT/pipeline_status.txt and $LOG_DIR" >&2
    exit "$overall"
fi

echo "Pipeline completed successfully: $RUN_ROOT"
