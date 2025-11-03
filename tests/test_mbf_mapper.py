from pathlib import Path

import numpy as np
import pytest

from mbf_mapper import map_mbf


# regression test to ensure that the results from the same input are the same
def test_map_mpf():
    dcm_series_fpath = Path("./data/MotionCorrectedPerfusionSeries/")
    tiff_masks_fpath = Path("./data/AIF_And_Myo_Masks.tiff")

    mbf_map_np = map_mbf(dcm_series_fpath, tiff_masks_fpath, is_vis_check=True)

    val_mean = np.nanmean(mbf_map_np)
    val_specific = mbf_map_np[105, 100]

    assert val_mean == pytest.approx(25.941189, abs=1e-5)
    assert val_specific == pytest.approx(23.781065, abs=1e-5)
