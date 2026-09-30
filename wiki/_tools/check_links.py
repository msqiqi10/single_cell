#!/usr/bin/env python3
"""Check a GitHub-wiki folder.

1. Broken [[...]] links: [[Page]] or [[Label|Page]] must resolve to Page.md.
2. Pages not reachable from Home.md or _Sidebar.md (following [[...]] links).
3. (with --paths) inline-code repo paths that do not exist relative to the repo
   root (the parent of the wiki folder). Wildcards and '...' are skipped.

Usage: python3 wiki/_tools/check_links.py [--paths]
Exit code 1 if any problem is found.
"""
import re
import sys
from pathlib import Path

WIKI = Path(__file__).resolve().parents[1]
ROOT = WIKI.parent
FENCE = re.compile(r"```.*?```", re.S)
INLINE = re.compile(r"`([^`\n]+)`")
LINK = re.compile(r"\[\[([^\]\n]+)\]\]")


def norm(t):
    return t.strip().replace(" ", "-")


def links_of(text):
    text = FENCE.sub("", text)
    text = INLINE.sub("", text)
    out = []
    for m in LINK.finditer(text):
        target = m.group(1).split("|")[-1]
        out.append(norm(target))
    return out


def main():
    pages = {p.stem: p for p in WIKI.glob("*.md") if p.stem != "README_PUBLISHING"}  # repo-only file, not a wiki page
    problems = 0
    graph = {}
    for name, p in pages.items():
        ls = links_of(p.read_text(encoding="utf-8"))
        graph[name] = ls
        for t in ls:
            if t not in pages:
                print(f"BROKEN  {name}.md -> [[{t}]]")
                problems += 1
    seen, stack = set(), [s for s in ("Home", "_Sidebar") if s in pages]
    seen.update(stack)
    while stack:
        cur = stack.pop()
        for t in graph.get(cur, []):
            if t in pages and t not in seen:
                seen.add(t)
                stack.append(t)
    for name in sorted(pages):
        if name not in seen and not name.startswith("_"):
            print(f"UNREACHABLE  {name}.md")
            problems += 1
    concept_missing = 0
    for name, p in pages.items():
        txt = p.read_text(encoding="utf-8")
        if "## 在本项目中" in txt and "## 相关概念" not in txt:
            print(f"MISSING 相关概念  {name}.md")
            concept_missing += 1
        if "## 相关概念" in txt and "## 在本项目中" not in txt:
            print(f"MISSING 在本项目中  {name}.md")
            concept_missing += 1
    problems += concept_missing
    if "--paths" in sys.argv:
        for name, p in sorted(pages.items()):
            txt = FENCE.sub("", p.read_text(encoding="utf-8"))
            for m in INLINE.finditer(txt):
                s = m.group(1).strip()
                if not re.match(r"^(/|[A-Za-z0-9_\-]+/)", s):
                    continue
                if not re.search(r"\.(md|csv|gz|json|py|pdf|pptx|html|txt|tsv|sh|npz|vtt|jpg|xlsx)$", s):
                    continue
                if "*" in s or "..." in s or " " in s.split("/")[-1] and False:
                    continue
                path = Path(s) if s.startswith("/") else ROOT / s
                if not path.exists():
                    print(f"PATH?  {name}.md -> {s}")
                    problems += 1
    print(f"pages: {len([n for n in pages if not n.startswith('_')])} (+ _Sidebar/_Footer), problems: {problems}")
    sys.exit(1 if problems else 0)


if __name__ == "__main__":
    main()
