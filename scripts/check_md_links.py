#!/usr/bin/env python3
"""Check relative links in all tracked-style markdown files of the repo.

Checks inline links [text](target) and <a href="target"> whose target is a
relative path (http/https/mailto/#anchor are skipped). Wiki [[...]] links are
covered by wiki/_tools/check_links.py. Optional '#fragment' and '?query' are
stripped. Fenced and inline code is ignored.

Absolute paths (server paths such as /home/zzz0054/...) are not checked.
Links listed in scripts/md_links_known_dangling.txt point at files that were
already absent from the 2026-09-29 download snapshot; they are reported as
"known" and do not fail the run. Exit 1 only on NEW broken links.

Usage: python3 scripts/check_md_links.py
"""
import re
import sys
from pathlib import Path
from urllib.parse import unquote

ROOT = Path(__file__).resolve().parents[1]
SKIP_DIRS = {".git", ".venv", "node_modules", "__pycache__"}
FENCE = re.compile(r"```.*?```|~~~.*?~~~", re.S)
INLINE = re.compile(r"`[^`\n]*`")
LINK = re.compile(r"(?<!\!)\[[^\]\n]*\]\(\s*<?([^)\s>]+)>?(?:\s+\"[^\"]*\")?\s*\)|!\[[^\]\n]*\]\(\s*<?([^)\s>]+)>?[^)]*\)|href=[\"']([^\"']+)[\"']")


def main():
    kf = ROOT / "scripts/md_links_known_dangling.txt"
    known = set(kf.read_text().splitlines()) if kf.exists() else set()
    bad = total = nknown = 0
    for md in sorted(ROOT.rglob("*.md")):
        if SKIP_DIRS & set(md.relative_to(ROOT).parts):
            continue
        text = INLINE.sub("", FENCE.sub("", md.read_text(encoding="utf-8", errors="replace")))
        for m in LINK.finditer(text):
            tgt = next(g for g in m.groups() if g)
            if re.match(r"^(https?:|mailto:|#|data:|tel:)", tgt) or tgt.startswith("[["):
                continue
            total += 1
            path = unquote(re.split(r"[#?]", tgt)[0])
            if not path:
                continue
            if path.startswith("/"):
                continue
            if not (md.parent / path).exists():
                if f"{md.relative_to(ROOT)} -> {tgt}" in known:
                    nknown += 1
                    continue
                bad += 1
                print(f"BROKEN {md.relative_to(ROOT)} -> {tgt}")
    print(f"links checked: {total}, known-dangling: {nknown}, NEW broken: {bad}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
