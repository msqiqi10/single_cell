#!/usr/bin/env bash
set -euo pipefail

ROOT="/home/zzz0054/bio3"
PY="$ROOT/.venv/bin/python"
SIG="$ROOT/notebooks/scripts/iNKT/analyze_inkt_gene_signatures.py"
SIG_TEST="$ROOT/notebooks/scripts/iNKT/test_analyze_inkt_gene_signatures.py"
DE="$ROOT/notebooks/scripts/iNKT/run_inkt_de_pathway_overlap.py"
TRAJ="$ROOT/notebooks/scripts/iNKT/analyze_inkt_trajectory_biomarkers.py"

cd "$ROOT"

echo "[1/4] Applying idempotent fixes..."

"$PY" - <<'PY'
from pathlib import Path

root = Path("/home/zzz0054/bio3")
sig = root / "notebooks/scripts/iNKT/analyze_inkt_gene_signatures.py"
test = root / "notebooks/scripts/iNKT/test_analyze_inkt_gene_signatures.py"
traj = root / "notebooks/scripts/iNKT/analyze_inkt_trajectory_biomarkers.py"
de = root / "notebooks/scripts/iNKT/run_inkt_de_pathway_overlap.py"

def swap(path, old, new):
    text = path.read_text()
    if new in text:
        print("already fixed:", path.name)
        return
    if old not in text:
        raise RuntimeError(f"Expected block not found in {path}")
    path.write_text(text.replace(old, new, 1))
    print("fixed:", path.name)

swap(
    test,
    "import importlib.util\nimport tempfile",
    "import importlib.util\nimport sys\nimport tempfile",
)
swap(
    test,
    "assert SPEC.loader is not None\nSPEC.loader.exec_module(MODULE)",
    "assert SPEC.loader is not None\nsys.modules[SPEC.name] = MODULE\nSPEC.loader.exec_module(MODULE)",
)

swap(
    sig,
    '''    match = re.match(r"(inkt1|inkt2|inkt17) \\(wang et al\\.,? 2022\\)", lower)
    if match:
''',
    '''    match = re.match(r"(inkt1|inkt2|inkt17) \\(wang et al\\.,? 2022\\)", lower)
    if match is None:
        match = re.fullmatch(r"wang\\s*2022\\s+(inkt1|inkt2|inkt17)", lower)
    if match:
''',
)
swap(sig, '    if lower.startswith("circulatory markers"):', '    if lower.startswith("circulatory"):')
swap(sig, '    if lower.startswith("residency markers"):', '    if lower.startswith("residency"):')

swap(
    sig,
    '''        union = sorted(
            set().union(*(set(canonical_sets[descriptor.set_id]) for descriptor in selected)),
            key=str.upper,
        )
''',
    '''        representative: dict[str, str] = {}
        for descriptor in selected:
            for gene in canonical_sets[descriptor.set_id]:
                representative.setdefault(gene.upper(), gene)
        union = [representative[key] for key in sorted(representative)]
''',
)

swap(
    sig,
    '''    id_columns = [
        column
        for column in summary.columns
        if column
        not in {
            "condition",
            "n_cells",
            *value_columns,
        }
    ]
''',
    '''    measurement_columns = set(value_columns) | {"std_score"}
    id_columns = [
        column
        for column in summary.columns
        if column not in {"condition", "n_cells", *measurement_columns}
    ]
''',
)

swap(
    sig,
    '        "input_h5ad_sha256": sha256_file(args.h5ad),',
    '        "input_h5ad_size_bytes": args.h5ad.stat().st_size,',
)

swap(
    traj,
    '''        actual = available.get(wanted.casefold())
        if actual is None:
            candidates = [name for name in book.sheet_names if name.casefold().startswith(wanted.casefold())]
            if len(candidates) != 1:
                raise KeyError(f"Could not uniquely resolve sheet {wanted!r}; found {candidates}")
            actual = candidates[0]
''',
    '''        actual = available.get(wanted.casefold())
        if actual is None:
            if "wang2022" in program:
                subtype = program.split("_", 1)[0]
                subtype_pattern = re.compile(re.escape(subtype) + r"(?!\\d)")
                candidates = [
                    name for name in book.sheet_names
                    if subtype_pattern.search(name.casefold()) and "wang" in name.casefold()
                ]
            else:
                candidates = [
                    name for name in book.sheet_names
                    if name.casefold().startswith(wanted.casefold())
                ]
            if len(candidates) != 1:
                raise KeyError(f"Could not uniquely resolve sheet {wanted!r}; found {candidates}")
            actual = candidates[0]
''',
)

swap(
    traj,
    '''        cell_scores.loc[finite_pt, "dpt_bin"] = pd.qcut(
            pt[finite_pt], q=args.n_bins, labels=False, duplicates="drop"
        ).astype("Int64").to_numpy()
''',
    '''        bin_values = pd.qcut(
            pd.Series(pt[finite_pt]), q=args.n_bins, labels=False, duplicates="drop"
        )
        cell_scores.loc[finite_pt, "dpt_bin"] = bin_values.astype("Int64").to_numpy()
''',
)

swap(
    de,
    '''    frames = [pd.read_csv(path, low_memory=False) for path in paths]
    return pd.concat(frames, ignore_index=True)
''',
    '''    frames = [pd.read_csv(path, low_memory=False) for path in paths]
    combined = pd.concat(frames, ignore_index=True)
    if "cluster_current" in combined:
        combined["cluster_current"] = combined["cluster_current"].map(
            lambda value: None if pd.isna(value) else str(int(float(value)))
        )
    return combined
''',
)
PY

echo "[2/4] Compiling and validating..."

"$PY" -m py_compile \
  "$SIG" \
  "$SIG_TEST" \
  "$DE" \
  "$TRAJ"

"$PY" "$SIG_TEST" -v
"$PY" "$SIG" --validate-only
"$PY" "$TRAJ" --validate-only
"$PY" "$DE" --help >/dev/null

echo "[3/4] Preparing parallel tmux jobs..."

STAMP="$(date -u +%Y%m%d_%H%M%S)"
SESSION="inkt_xlsx_${STAMP}"
ANALYSIS_ROOT="$ROOT/output/iNKT_extended_runs/$STAMP"
TMUX_DIR="$ANALYSIS_ROOT/tmux"

mkdir -p "$TMUX_DIR/jobs" "$TMUX_DIR/logs" "$TMUX_DIR/status"

make_job() {
    local name="$1"
    local threads="$2"
    shift 2

    local job="$TMUX_DIR/jobs/${name}.sh"
    local log="$TMUX_DIR/logs/${name}.log"
    local status="$TMUX_DIR/status/${name}.exit"

    {
        printf '#!/usr/bin/env bash\n'
        printf 'set -u -o pipefail\n'
        printf 'cd %q\n' "$ROOT"
        printf 'export CUDA_VISIBLE_DEVICES=""\n'
        printf 'export OMP_NUM_THREADS=%q\n' "$threads"
        printf 'export OPENBLAS_NUM_THREADS=%q\n' "$threads"
        printf 'export MKL_NUM_THREADS=%q\n' "$threads"
        printf 'export NUMEXPR_NUM_THREADS=%q\n' "$threads"
        printf '{ date -u +"started_utc=%%FT%%TZ"; '
        printf '%q ' "$@"
        printf '; } > %q 2>&1\n' "$log"
        printf 'rc=$?\n'
        printf 'printf "%%s\\n" "$rc" > %q\n' "$status"
        printf 'printf "finished_utc=%%s rc=%%s\\n" "$(date -u +%%FT%%TZ)" "$rc"\n'
        printf 'exit "$rc"\n'
    } > "$job"
}

make_job gene_sets 1 \
  "$PY" "$SIG" \
  --task gene-sets --threads 1 \
  --out-dir "$ANALYSIS_ROOT/signatures"

make_job signature_scores 3 \
  "$PY" "$SIG" \
  --task signature-scores --threads 3 \
  --out-dir "$ANALYSIS_ROOT/signatures"

make_job marker_validation 2 \
  "$PY" "$SIG" \
  --task marker-validation --threads 2 \
  --out-dir "$ANALYSIS_ROOT/signatures"

make_job legacy 1 \
  "$PY" "$DE" \
  --task legacy \
  --outdir "$ANALYSIS_ROOT/de_pathway"

make_job de0 1 \
  "$PY" "$DE" \
  --task current-de --n-jobs 7 \
  --shard-count 2 --shard-index 0 \
  --outdir "$ANALYSIS_ROOT/de_pathway"

make_job de1 1 \
  "$PY" "$DE" \
  --task current-de --n-jobs 7 \
  --shard-count 2 --shard-index 1 \
  --outdir "$ANALYSIS_ROOT/de_pathway"

make_job trajectory 3 \
  "$PY" "$TRAJ" \
  --output-dir "$ANALYSIS_ROOT/trajectory"

cat > "$TMUX_DIR/jobs/finalize.sh" <<'FINALIZE'
#!/usr/bin/env bash
set -u -o pipefail

DE_SCRIPT="$ROOT/notebooks/scripts/iNKT/run_inkt_de_pathway_overlap.py"
PYTHON="$ROOT/.venv/bin/python"
STATUS="$RUN_DIR/tmux/status"
LOGS="$RUN_DIR/tmux/logs"

deps=(legacy de0 de1)

echo "Waiting for legacy and both DE shards..."
while true; do
    ready=1
    for name in "${deps[@]}"; do
        [[ -f "$STATUS/${name}.exit" ]] || ready=0
    done
    [[ "$ready" -eq 1 ]] && break
    date -u +"waiting_utc=%FT%TZ"
    sleep 15
done

dependency_ok=1
for name in "${deps[@]}"; do
    code="$(<"$STATUS/${name}.exit")"
    [[ "$code" == "0" ]] || dependency_ok=0
done

if [[ "$dependency_ok" -eq 1 ]]; then
    echo "Dependencies passed; running overlap/pathway integration..."
    env \
      CUDA_VISIBLE_DEVICES="" \
      OMP_NUM_THREADS=1 \
      OPENBLAS_NUM_THREADS=1 \
      MKL_NUM_THREADS=1 \
      "$PYTHON" "$DE_SCRIPT" \
      --task overlap-pathway \
      --outdir "$RUN_DIR/de_pathway" \
      >"$LOGS/overlap_pathway.log" 2>&1
    overlap_rc=$?
else
    echo "Skipping overlap/pathway because a dependency failed."
    overlap_rc=90
fi

printf '%s\n' "$overlap_rc" > "$STATUS/overlap_pathway.exit"

all_jobs=(
  gene_sets
  signature_scores
  marker_validation
  legacy
  de0
  de1
  trajectory
  overlap_pathway
)

echo "Waiting for all jobs..."
while true; do
    ready=1
    for name in "${all_jobs[@]}"; do
        [[ -f "$STATUS/${name}.exit" ]] || ready=0
    done
    [[ "$ready" -eq 1 ]] && break
    sleep 15
done

overall=0
{
    date -u +"completed_utc=%FT%TZ"
    for name in "${all_jobs[@]}"; do
        code="$(<"$STATUS/${name}.exit")"
        printf '%s=%s\n' "$name" "$code"
        [[ "$code" == "0" ]] || overall=1
    done
    printf 'overall=%s\n' "$overall"
} | tee "$RUN_DIR/pipeline_status.txt"

printf '%s\n' "$overall" > "$STATUS/pipeline.exit"
exit "$overall"
FINALIZE

echo "[4/4] Starting tmux session $SESSION ..."

tmux new-session -d -s "$SESSION" -n bootstrap "sleep 10"
tmux set-option -t "$SESSION" remain-on-exit on

for name in \
    gene_sets \
    signature_scores \
    marker_validation \
    legacy \
    de0 \
    de1 \
    trajectory
do
    tmux new-window \
      -t "$SESSION" \
      -n "$name" \
      "bash '$TMUX_DIR/jobs/${name}.sh'"
done

tmux new-window \
  -t "$SESSION" \
  -n finalize \
  "env ROOT='$ROOT' RUN_DIR='$ANALYSIS_ROOT' bash '$TMUX_DIR/jobs/finalize.sh' >'$TMUX_DIR/logs/finalize.log' 2>&1"

tmux kill-window -t "$SESSION:bootstrap"

printf '%s\n' "$SESSION" > "$ROOT/.inkt_tmux_session"
printf '%s\n' "$ANALYSIS_ROOT" > "$ROOT/.inkt_tmux_run_dir"

echo
echo "Started successfully."
echo "Session: $SESSION"
echo "Outputs: $ANALYSIS_ROOT"
echo
echo "Attach: tmux attach -t $SESSION"
echo "Detach: Ctrl-b then d"
echo "Progress: tail -f $TMUX_DIR/logs/finalize.log"
