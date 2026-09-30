from __future__ import annotations

from pathlib import Path

import nbformat
from nbclient import NotebookClient


ROOT = Path(__file__).resolve().parents[3]
NOTEBOOK = ROOT / "notebooks/iNKT/scanpy_iNKT_preprocess_plotting_trajectory.ipynb"

with NOTEBOOK.open() as handle:
    notebook = nbformat.read(handle, as_version=4)

client = NotebookClient(
    notebook,
    timeout=7200,
    kernel_name="python3",
    allow_errors=False,
    resources={"metadata": {"path": str(ROOT)}},
)
client.execute()

with NOTEBOOK.open("w") as handle:
    nbformat.write(notebook, handle)

print(NOTEBOOK)
