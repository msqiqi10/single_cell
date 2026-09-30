from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

import pandas as pd


# Keep ingestion order stable so a palette-only change cannot perturb row order
# and, in turn, the computed embedding.
SAMPLE_LOAD_ORDER = (
    "Ctrl_BM",
    "Ctrl_Spleen",
    "Ctrl_Thymus",
    "T2_BM",
    "T2_Spleen",
    "T2_Thymus",
)

# The legacy PPT groups samples by condition and orders tissues Thymus, BM,
# Spleen. Tissue controls hue; condition selects the light/dark hue variant.
SAMPLE_DISPLAY_ORDER = (
    "Ctrl_Thymus",
    "Ctrl_BM",
    "Ctrl_Spleen",
    "T2_Thymus",
    "T2_BM",
    "T2_Spleen",
)
SAMPLE_PALETTE = {
    "Ctrl_Thymus": "#0000FF",
    "Ctrl_BM": "#FFFF00",
    "Ctrl_Spleen": "#FF0000",
    "T2_Thymus": "#00008B",
    "T2_BM": "#FFA500",
    "T2_Spleen": "#8B0000",
}

CONDITION_ORDER = ("Ctrl", "T2")
# The PPT has no global two-color condition palette. These derived neutral
# colors preserve its light/dark condition logic without borrowing a tissue hue.
CONDITION_PALETTE = {
    "Ctrl": "#A6A6A6",
    "T2": "#4D4D4D",
}

TISSUE_MAP = {
    "BM": "bone_marrow",
    "Spleen": "spleen",
    "Thymus": "thymus",
}
TISSUE_ORDER = ("thymus", "bone_marrow", "spleen")
TISSUE_PALETTE = {
    "thymus": "#0000FF",
    "bone_marrow": "#FFFF00",
    "spleen": "#FF0000",
}

PALETTE_SPECS: dict[str, tuple[tuple[str, ...], dict[str, str]]] = {
    "sample": (SAMPLE_DISPLAY_ORDER, SAMPLE_PALETTE),
    "condition": (CONDITION_ORDER, CONDITION_PALETTE),
    "tissue": (TISSUE_ORDER, TISSUE_PALETTE),
}


def _validate_palette_definition(
    key: str,
    categories: Sequence[str],
    palette: Mapping[str, str],
) -> None:
    category_set = set(categories)
    palette_set = set(palette)
    if len(category_set) != len(categories):
        raise ValueError(f"Duplicate categories in {key} order: {categories}")
    if palette_set != category_set:
        missing = sorted(category_set - palette_set)
        extra = sorted(palette_set - category_set)
        raise ValueError(f"Invalid {key} palette; missing={missing}, extra={extra}")


def set_categorical_palette(
    adata: Any,
    key: str,
    categories: Sequence[str],
    palette: Mapping[str, str],
) -> None:
    """Validate metadata, set display order, and install a Scanpy palette."""
    _validate_palette_definition(key, categories, palette)
    if key not in adata.obs:
        raise KeyError(f"Missing obs column: {key}")

    values = adata.obs[key]
    if values.isna().any():
        raise ValueError(f"obs[{key!r}] contains missing values before categorization")
    string_values = [str(value) for value in values]
    unknown = sorted(set(string_values) - set(categories))
    if unknown:
        raise ValueError(f"Unknown {key} values: {unknown}")

    adata.obs[key] = pd.Categorical(string_values, categories=list(categories))
    if adata.obs[key].isna().any():
        raise ValueError(f"obs[{key!r}] contains missing values after categorization")
    adata.uns[f"{key}_colors"] = [palette[category] for category in categories]


def configure_inkt_palettes(adata: Any) -> None:
    for key, (categories, palette) in PALETTE_SPECS.items():
        set_categorical_palette(adata, key, categories, palette)


def palette_state(adata: Any) -> dict[str, dict[str, list[str]]]:
    state: dict[str, dict[str, list[str]]] = {}
    for key in PALETTE_SPECS:
        series = adata.obs.get(key)
        if series is not None and isinstance(series.dtype, pd.CategoricalDtype):
            categories = [str(value) for value in series.cat.categories]
        else:
            categories = []
        raw_colors = adata.uns.get(f"{key}_colors", [])
        if hasattr(raw_colors, "tolist"):
            raw_colors = raw_colors.tolist()
        state[key] = {
            "categories": categories,
            "colors": [str(value).upper() for value in raw_colors],
        }
    return state


def palette_validation_failures(adata: Any) -> list[str]:
    failures: list[str] = []
    for key, (expected_categories, palette) in PALETTE_SPECS.items():
        if key not in adata.obs:
            failures.append(f"obs column is missing: {key}")
            continue
        series = adata.obs[key]
        if not isinstance(series.dtype, pd.CategoricalDtype):
            failures.append(f"obs column is not categorical: {key}")
        else:
            actual_categories = tuple(str(value) for value in series.cat.categories)
            if actual_categories != expected_categories:
                failures.append(
                    f"{key} categories {actual_categories} != {expected_categories}"
                )
        if series.isna().any():
            failures.append(f"obs column contains missing values: {key}")

        uns_key = f"{key}_colors"
        if uns_key not in adata.uns:
            failures.append(f"uns entry is missing: {uns_key}")
            continue
        raw_colors = adata.uns[uns_key]
        if hasattr(raw_colors, "tolist"):
            raw_colors = raw_colors.tolist()
        actual_colors = tuple(str(value).upper() for value in raw_colors)
        expected_colors = tuple(palette[category].upper() for category in expected_categories)
        if actual_colors != expected_colors:
            failures.append(f"{uns_key} {actual_colors} != {expected_colors}")
    return failures


for _key, (_categories, _palette) in PALETTE_SPECS.items():
    _validate_palette_definition(_key, _categories, _palette)
