from pathlib import Path
import re
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'2026-09-30'; RES=OUT/'results'
P19=ROOT/'2026-09-19'; P25=ROOT/'2026-09-25'
KEGG=Path('/Users/zeruzhang/Documents/Codex/2026-09-24/users-zeruzhang-downloads-zoom-0923/work/remote37/probe/output/iNKT_reproduction_deck/20260830_C5_paper_Fig3DEF_followup/tables/20260830_KEGG_2019_Mouse.gmt')
MINILM=Path('/Users/zeruzhang/Documents/Codex/2026-09-24/users-zeruzhang-downloads-zoom-0923/local_inkt/models/all-MiniLM-L6-v2')
VENV_PY='/Users/zeruzhang/Documents/Codex/2026-09-24/users-zeruzhang-downloads-zoom-0923/local_inkt/.venv/bin/python'
# Version (b) exclusion. Hsp(?!g) keeps Hspg2 (perlecan, not a heat-shock protein); Rps6k* kinases (MAPK-relevant) are NOT matched.
EXCL_REGEX=r'^(Rp[ls]\d|Rplp\d|Rpsa$|Mrp[ls]\d|mt-|Hsp(?!g)|Dnaj)'
EXCL=re.compile(EXCL_REGEX)
ANCHOR_KEGG=['MAPK signaling pathway','IL-17 signaling pathway','T cell receptor signaling pathway','Th17 cell differentiation','TNF signaling pathway']
ANCHOR_GO=['GO:0035976','MAPK cascade','ERK1 and ERK2 cascade','p38MAPK cascade','JNK cascade','T-helper 17 type immune response','interleukin-17 production']
