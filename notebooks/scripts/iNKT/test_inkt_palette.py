from __future__ import annotations

import unittest

import anndata as ad
import numpy as np
import pandas as pd

from inkt_palette import (
    CONDITION_ORDER,
    CONDITION_PALETTE,
    SAMPLE_DISPLAY_ORDER,
    SAMPLE_LOAD_ORDER,
    SAMPLE_PALETTE,
    TISSUE_ORDER,
    TISSUE_PALETTE,
    configure_inkt_palettes,
    palette_state,
    palette_validation_failures,
)


def example_adata() -> ad.AnnData:
    obs = pd.DataFrame(
        {
            "sample": list(SAMPLE_LOAD_ORDER),
            "condition": [value.split("_", 1)[0] for value in SAMPLE_LOAD_ORDER],
            "tissue": [
                {
                    "BM": "bone_marrow",
                    "Spleen": "spleen",
                    "Thymus": "thymus",
                }[value.split("_", 1)[1]]
                for value in SAMPLE_LOAD_ORDER
            ],
        },
        index=[f"cell_{idx}" for idx in range(len(SAMPLE_LOAD_ORDER))],
    )
    return ad.AnnData(X=np.zeros((len(obs), 1)), obs=obs)


class InktPaletteTests(unittest.TestCase):
    def test_canonical_palette_matches_ppt(self) -> None:
        self.assertEqual(
            SAMPLE_PALETTE,
            {
                "Ctrl_Thymus": "#0000FF",
                "Ctrl_BM": "#FFFF00",
                "Ctrl_Spleen": "#FF0000",
                "T2_Thymus": "#00008B",
                "T2_BM": "#FFA500",
                "T2_Spleen": "#8B0000",
            },
        )
        self.assertEqual(tuple(TISSUE_PALETTE), TISSUE_ORDER)
        self.assertEqual(tuple(CONDITION_PALETTE), CONDITION_ORDER)

    def test_configure_preserves_rows_and_sets_exact_ordered_palettes(self) -> None:
        adata = example_adata()
        original_rows = adata.obs_names.tolist()
        original_values = adata.obs["sample"].tolist()

        configure_inkt_palettes(adata)

        self.assertEqual(adata.obs_names.tolist(), original_rows)
        self.assertEqual(adata.obs["sample"].astype(str).tolist(), original_values)
        self.assertEqual(tuple(adata.obs["sample"].cat.categories), SAMPLE_DISPLAY_ORDER)
        self.assertEqual(tuple(adata.obs["condition"].cat.categories), CONDITION_ORDER)
        self.assertEqual(tuple(adata.obs["tissue"].cat.categories), TISSUE_ORDER)
        self.assertEqual(palette_validation_failures(adata), [])
        self.assertEqual(
            palette_state(adata)["sample"]["colors"],
            [SAMPLE_PALETTE[category] for category in SAMPLE_DISPLAY_ORDER],
        )

    def test_unknown_metadata_fails_fast(self) -> None:
        adata = example_adata()
        adata.obs.loc[adata.obs.index[0], "sample"] = "Ctrl_Liver"
        with self.assertRaisesRegex(ValueError, "Unknown sample values"):
            configure_inkt_palettes(adata)

    def test_palette_validation_detects_color_drift(self) -> None:
        adata = example_adata()
        configure_inkt_palettes(adata)
        adata.uns["tissue_colors"] = ["#000000"] * len(TISSUE_ORDER)
        failures = palette_validation_failures(adata)
        self.assertTrue(any("tissue_colors" in failure for failure in failures))


if __name__ == "__main__":
    unittest.main()
