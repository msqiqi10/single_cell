from __future__ import annotations

import argparse
from pathlib import Path

import nbformat
from nbclient import NotebookClient


ROOT = Path(__file__).resolve().parents[3]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Execute an iNKT notebook in the repository environment.")
    parser.add_argument("--notebook", type=Path, required=True)
    parser.add_argument("--timeout", type=int, default=7200)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    notebook_path = args.notebook.resolve()
    if not notebook_path.exists():
        raise FileNotFoundError(notebook_path)
    if args.timeout < 1:
        raise ValueError("--timeout must be >= 1")

    notebook = nbformat.read(notebook_path, as_version=4)
    client = NotebookClient(
        notebook,
        timeout=args.timeout,
        kernel_name="python3",
        allow_errors=False,
        resources={"metadata": {"path": str(ROOT)}},
    )
    client.execute()
    nbformat.write(notebook, notebook_path)
    print(notebook_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
