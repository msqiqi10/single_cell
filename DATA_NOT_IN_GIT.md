# Data not in Git

The files below are kept on disk but excluded by `.gitignore` (`*.h5ad`, `*.zip`).
Git LFS is deliberately not set up. GitHub rejects files over 100 MB and warns above 50 MB.
Three of these files exceed 50 MB (marked **>50 MB**); the other four are small or mid-sized
but are ignored by the same `*.h5ad` rule. Total: 430,334,716 bytes (about 410 MiB).

Original location of everything: server `/home/zzz0054/bio3` (same relative paths).
SHA-256 values were computed locally on 2026-09-30 and match `source_manifest.json` from the 2026-09-29 download.

| Path | Size (bytes) | SHA-256 | What it is | Regenerate / obtain |
|---|---:|---|---|---|
| `output/iNKT_meeting_followup_20260905/objects/scored_base.h5ad` **>50 MB** | 213,444,986 | `f2c7e3c019044a48e8a08c08913e8617c77e6dc3abda432bf99100fd21cdb453` | AnnData of all 15,532 cells x 10,670 genes with module scores, UMAP and refined clusters | Server copy; rebuild with `scripts/runners/run_inkt_meeting_followup_20260905.sh scores` (runs `notebooks/scripts/iNKT/run_inkt_meeting_followup.py`) |
| `output/iNKT_meeting_followup_20260905/objects/bone_marrow_stable.h5ad` **>50 MB** | 84,211,048 | `57110aab150449ff1041a18f32d7de6b022c02f7106801a8301bfa789121ade2` | Bone-marrow subset with stable re-clustering (sensitivity analysis) | Same runner (stable-clustering stage) |
| `output/iNKT_meeting_followup_20260905/objects/spleen_stable.h5ad` **>50 MB** | 90,886,903 | `183ed066c26534f19928ca585ab010ca37abb113e47092245f7b016b266bd292` | Spleen subset with stable re-clustering | Same runner |
| `output/iNKT_meeting_followup_20260905/objects/thymus_stable.h5ad` | 39,354,471 | `e77b527ecb8fd7b245f1d4f9f46408d6e6526bcd0be362f641c841ce1dd00fd6` | Thymus subset with stable re-clustering | Same runner |
| `iNKT_by_date/2026-09-19/results/objects/bone_marrow_local_trajectory.h5ad` | 938,554 | `3600cbf33d3ea262af6132418fa4c7a74f301a547efb547359e7b90189af1414` | Bone-marrow local DPT trajectory object | `iNKT_by_date/2026-09-19/code/analyze.py` (needs scored_base.h5ad) |
| `iNKT_by_date/2026-09-19/results/objects/spleen_local_trajectory.h5ad` | 1,062,131 | `3ce0c519316bb961c7d821d8b132a1e01d1f74c6c44841e71c3bb7560e560c78` | Spleen local DPT trajectory object | same |
| `iNKT_by_date/2026-09-19/results/objects/thymus_local_trajectory.h5ad` | 436,623 | `dd933d59b9c8f82e4015b90e508cfbde243beb9c1027492db7a1b31aa131541b` | Thymus local DPT trajectory object | same |

Other things mentioned in the old README but never present in this snapshot: `input.zip` (~7 GB raw data archive),
`CML_NK_scRNA_TKI/` (shallow clone of `https://github.com/ai-pharm-AU/CML_NK_scRNA_TKI`), and `.venv/`
(recreate with `uv sync`). They are covered by `.gitignore` where applicable.

To restore: the 7 `.h5ad` files are hosted on the **private** Hugging Face dataset
[`si3g/inkt-scrna-data`](https://huggingface.co/datasets/si3g/inkt-scrna-data) under the same relative paths
(manifest with SHA-256: `data_manifest.json`). The repo is private, so access must be granted by the owner first.

```
hf auth login
python scripts/fetch_data.py          # download missing files, verify SHA-256
python scripts/fetch_data.py --check  # verify local files only
```

The script skips files that already have the right hash and never overwrites a local file whose hash differs.
Alternatively copy from the server (`scp` / `rsync` from `/home/zzz0054/bio3`) and verify with `shasum -a 256`.
