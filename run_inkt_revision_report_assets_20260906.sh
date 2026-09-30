#!/usr/bin/env bash
set -euo pipefail
cd /home/zzz0054/bio3
out=output/iNKT_revision_responses_20260906
trap 'echo $? > "$out/logs/assets.exit"' EXIT
/home/zzz0054/miniforge3/bin/uv pip install --python .venv/bin/python --target /tmp/inkt_revision_docx_tools python-docx
curl --fail --location --connect-timeout 20 --max-time 120 --output "$out/assets/NotoSansCJKsc-Regular.otf" https://raw.githubusercontent.com/notofonts/noto-cjk/main/Sans/OTF/SimplifiedChinese/NotoSansCJKsc-Regular.otf
curl --fail --location --connect-timeout 20 --max-time 60 --output "$out/assets/NotoSansCJK-LICENSE.txt" https://raw.githubusercontent.com/notofonts/noto-cjk/main/Sans/LICENSE
