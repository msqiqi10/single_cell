#!/usr/bin/env python3
"""Build compact visual-QA contact sheets for the iNKT reproduction deck."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


THUMBNAIL = (640, 360)
GRID = (2, 4)
MARGIN = 18
LABEL_HEIGHT = 24


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--deck-dir", type=Path, required=True)
    args = parser.parse_args()
    deck_dir = args.deck_dir.resolve()
    manifest = json.loads((deck_dir / "deck_manifest.json").read_text(encoding="utf-8"))
    slides = manifest["slides"]
    page_capacity = GRID[0] * GRID[1]
    sheet_width = MARGIN + GRID[0] * (THUMBNAIL[0] + MARGIN)
    sheet_height = MARGIN + GRID[1] * (THUMBNAIL[1] + LABEL_HEIGHT + MARGIN)
    font = ImageFont.load_default()

    for batch_start in range(0, len(slides), page_capacity):
        sheet = Image.new("RGB", (sheet_width, sheet_height), "white")
        draw = ImageDraw.Draw(sheet)
        for slot, record in enumerate(slides[batch_start : batch_start + page_capacity]):
            row, column = divmod(slot, GRID[0])
            x = MARGIN + column * (THUMBNAIL[0] + MARGIN)
            y = MARGIN + row * (THUMBNAIL[1] + LABEL_HEIGHT + MARGIN)
            slide_path = deck_dir / "slides" / f"slide_{int(record['slide']):02d}.png"
            with Image.open(slide_path) as source:
                image = source.convert("RGB")
                image.thumbnail(THUMBNAIL, Image.Resampling.LANCZOS)
            paste_x = x + (THUMBNAIL[0] - image.width) // 2
            paste_y = y + (THUMBNAIL[1] - image.height) // 2
            sheet.paste(image, (paste_x, paste_y))
            draw.rectangle(
                (x, y, x + THUMBNAIL[0], y + THUMBNAIL[1]),
                outline="#CBD5E1",
                width=1,
            )
            label = f"Slide {int(record['slide']):02d}: {record['title']}"
            draw.text((x + 2, y + THUMBNAIL[1] + 5), label[:92], fill="#102A43", font=font)
        sheet_index = batch_start // page_capacity + 1
        output = deck_dir / f"qa_contact_{sheet_index}.jpg"
        sheet.save(output, quality=92, optimize=True)

    print(f"contact_sheets={(len(slides) + page_capacity - 1) // page_capacity}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
