import os

import cv2
import numpy as np
import SimpleITK as sitk
import pydicom


def get_perf_ts(dcm_series_path):
    """Reads timestamps from DICOM series with perfusion."""
    filenames = os.listdir(dcm_series_path)

    times = []
    for filename in filenames:
        ds = pydicom.dcmread(dcm_series_path / filename, stop_before_pixels=True)
        ts_str = getattr(ds, "AcquisitionTime", None) or getattr(ds, "ContentTime", None)
        if ts_str and len(ts_str) >= 6:
            # convert to seconds
            h, m, s = int(ts_str[ : 2]), int(ts_str[2 : 4]), float(ts_str[4 : ])
            times.append(h * 3600 + m * 60 + s)

    times = sorted(times)
    # print(times)
    # dt_s = np.mean(np.diff(times)) if len(times) > 1 else None
    print("time diffs, sec: ", np.diff(times))
    # print(f"Estimated frame spacing: {dt_s:.3f} s")

    times_np = np.asarray(times)
    times_np = times_np - times_np[0]

    return times_np


def read_perf_dcm(dcm_series_path):
    """Reads dicom series with perfusion."""

    # volume
    reader = sitk.ImageSeriesReader()
    files = reader.GetGDCMSeriesFileNames(dcm_series_path)
    reader.SetFileNames(files)
    img = reader.Execute()

    arr = sitk.GetArrayFromImage(img)
    vol_np = arr.astype(np.float32)

    return vol_np


def read_roi_masks(tiff_masks_path):
    """Reads RoI masks of left ventricle and myocardium."""

    # load all pages from the TIFF
    success, images = cv2.imreadmulti(tiff_masks_path, [], cv2.IMREAD_UNCHANGED)

    if not success:
        raise ValueError("Could not read TIFF file")

    print(f"Number of masks extracted from {tiff_masks_path}: {len(images)}")
    
    mask_ventr = images[0]
    mask_myo = images[1]

    print(f"mask_myo.shape: {mask_myo.shape}")

    return mask_ventr, mask_myo