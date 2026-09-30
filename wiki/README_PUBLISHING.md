# Publishing this wiki to GitHub

GitHub wikis live in a separate repository, `<repo>.wiki.git`, not in the main repo.
This folder is written in GitHub-wiki layout (`Home.md`, `_Sidebar.md`, `_Footer.md`, one page per `.md`),
so it can be pushed there as-is. `_tools/` and this file are not wiki pages; delete them from the wiki copy
(or leave them, they are harmless but will show as pages).

Prerequisite: on GitHub, open the repo's Wiki tab and create any first page once, so that `<repo>.wiki.git` exists.

## Option A: copy the files (simplest)

```bash
git clone git@github.com:<user>/<repo>.wiki.git /tmp/wiki
rsync -a --exclude '_tools' --exclude 'README_PUBLISHING.md' wiki/ /tmp/wiki/
cd /tmp/wiki && git add -A && git commit -m "Sync wiki from main repo" && git push
```

## Option B: git subtree split (keeps the history of wiki/ only)

```bash
git subtree split --prefix=wiki -b wiki-only
git push git@github.com:<user>/<repo>.wiki.git wiki-only:master   # wiki default branch is usually master
```

Subsequent updates: repeat the split (it is incremental) and push, using `--force` only if the wiki
repo was also edited in the browser.

## Before publishing

- Run `python3 wiki/_tools/check_links.py --paths` from the repo root. It checks `[[wiki links]]` and that
  repo paths quoted in backticks (for example `iNKT_by_date/...`, `docs/audits/...`) still exist.
- Those inline paths are repo-relative and will appear as plain code text on the wiki, not clickable links.
  If you want clickable links, rewrite them to full `https://github.com/<user>/<repo>/blob/main/<path>` URLs.
- The wiki is the single source of truth only for its text; the numbers come from the repo files it cites.
