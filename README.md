# bio3 Scanpy Environment

This directory is configured as a `uv` project for Scanpy/scVI single-cell analysis.

## iNKT results by date

Start with the [dated iNKT directory](iNKT_by_date/README.md) for each stage's code,
results, presentations, and notes. The directory uses relative symbolic links to
the existing files; code links show the current source, not historical snapshots.

- [2026-09-15: reference-aligned presentation](iNKT_by_date/2026-09-15/README.md)
- [2026-09-06: R01–R13 report and complete delivery bundle](iNKT_by_date/2026-09-06/README.md)
- [2026-09-05: meeting follow-up and discovery results](iNKT_by_date/2026-09-05/README.md)
- [Experiment history and supervisor requirements](docs/inkt_history_and_supervisor_requirements_20260915.md)

## Layout

- `input.zip`: existing data archive. It has not been extracted.
- `CML_NK_scRNA_TKI/`: shallow clone of `https://github.com/ai-pharm-AU/CML_NK_scRNA_TKI`.
- `.venv/`: local virtual environment created by `uv sync`.

## Environment

Create or refresh the environment:

```bash
uv sync
```

Activate it:

```bash
source .venv/bin/activate
```

Smoke test:

```bash
uv run python -c "import scanpy as sc; import scvi; print(sc.__version__); print(scvi.__version__)"
```

The environment pins Python 3.12 and configures PyTorch from CPU wheels for Linux so setup does not target CUDA GPUs.

## Repo Notes

The cloned CML workflow code is in:

```text
CML_NK_scRNA_TKI/capstonebio-main/python/
```

Those scripts currently contain hard-coded absolute paths from the original author's machine. Update paths before running the workflows against local data.

If you decide to extract the data archive, check available disk first because `input.zip` is about 7 GB.
