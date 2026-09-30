#!/usr/bin/env python3
"""Fetch large h5ad files (not in git) from the private HF dataset si3g/inkt-scrna-data.

Usage: python scripts/fetch_data.py [--check]
Requires `hf auth login` (repo is private) and huggingface_hub for downloads.
Never overwrites a local file whose hash differs; it reports it instead.
"""
import argparse, hashlib, json, shutil, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true", help="only verify local files, no download")
    a = ap.parse_args()
    m = json.load(open(ROOT / "data_manifest.json"))
    bad = 0
    for e in m["files"]:
        dest = ROOT / e["path"]
        if dest.exists():
            if sha256(dest) == e["sha256"]:
                print(f"OK        {e['path']}")
            else:
                print(f"MISMATCH  {e['path']} (local file differs; not overwritten)")
                bad += 1
            continue
        if a.check:
            print(f"MISSING   {e['path']}")
            bad += 1
            continue
        from huggingface_hub import hf_hub_download
        print(f"DOWNLOAD  {e['path']}")
        tmp = Path(hf_hub_download(m["repo_id"], e["path"], repo_type=m["repo_type"]))
        if sha256(tmp) != e["sha256"]:
            print(f"BADHASH   {e['path']} (downloaded file hash mismatch; not installed)")
            bad += 1
            continue
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(tmp, dest)
        print(f"OK        {e['path']}")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
