#!/usr/bin/env bash
set -euo pipefail
cd /home/zzz0054/bio3
out=output/iNKT_revision_responses_20260906
export CUDA_VISIBLE_DEVICES=""
trap 'echo $? > "$out/logs/render.exit"' EXIT
/home/zzz0054/miniforge3/bin/uv pip install --python .venv/bin/python --target /tmp/inkt_revision_playwright_tools playwright
PYTHONPATH=/tmp/inkt_revision_playwright_tools .venv/bin/python -u - <<'PY'
from pathlib import Path
from playwright.sync_api import sync_playwright
p=Path('/home/zzz0054/bio3/output/iNKT_revision_responses_20260906')
with sync_playwright() as w:
 b=w.chromium.launch(executable_path='/home/zzz0054/.cache/ms-playwright/chromium-1217/chrome-linux64/chrome',headless=True,args=['--no-sandbox','--disable-gpu','--disable-dev-shm-usage'])
 page=b.new_page(viewport={'width':1250,'height':950});print('loading document',flush=True)
 page.goto((p/'iNKT_revision_responses_20260906.html').as_uri(),wait_until='load',timeout=45000)
 page.evaluate('document.fonts.ready');print('font and content loaded',flush=True)
 print(page.evaluate('({images:document.images.length,loaded:[...document.images].filter(i=>i.complete&&i.naturalWidth>0).length,fontStatus:document.fonts.status})'),flush=True)
 page.pdf(path=str(p/'iNKT_revision_responses_20260906.pdf'),format='A4',print_background=True,prefer_css_page_size=True)
 page.screenshot(path=str(p/'assets/report_browser_top.png'))
 b.close()
print('PDF complete',flush=True)
PY
