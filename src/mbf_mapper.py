import argparse
import os
from pathlib import Path

import cv2
import numpy as np
from numpy.fft import rfft, irfft
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
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

    print(f"Number of masks: {len(images)}")
    print(f"Shape of first mask: {images[0].shape}")

    mask_ventr = images[0]
    mask_myo = images[1]

    return mask_ventr, mask_myo


def get_aif_curve(vol, mask_bpool):
    """Returns discreet AIF curve from cardiac perfusion scan and corresponding blood pool (e.g. left ventricle) mask"""
    
    vol_bpool = np.where(mask_bpool, vol, 0)

    ts_num, _, _ = vol_bpool.shape

    # Robust mean over non-zero blood-pool pixels each frame
    aif_t = []
    for ts in range(ts_num):
        vals = vol_bpool[ts][vol_bpool[ts] > 0]
        aif_t.append(np.median(vals) if vals.size else 0.0)
    aif_t = np.asarray(aif_t, dtype=np.float64)

    return aif_t


def next_pow2(x):
    return 1 << (int(x - 1).bit_length())

def wiener_deconv_all(aif, Y_tN, dt, alpha=1e-2):
    """
    Vectorized Wiener deconvolution for many curves at once.
    aif: (T,)
    Y_tN: (T, N) tissue curves (columns = pixels)
    returns h_tN: (T, N)
    """
    T, N = Y_tN.shape
    n = next_pow2(2*T)  # zero-pad to avoid circular wrap-around

    a_pad = np.pad(aif, (0, n - T))
    A = rfft(a_pad)                              # (F,)
    denom = (np.abs(A)**2 + alpha)               # (F,)

    Y_pad = np.pad(Y_tN, ((0, n - T), (0, 0)))   # (n, N)
    YF = rfft(Y_pad, axis=0)                     # (F, N)

    HF = (YF * np.conj(A)[:, None]) / denom[:, None]   # (F, N)
    h_tN = irfft(HF, n=n, axis=0)[:T, :] / dt          # back to time, scale by dt
    h_tN = np.clip(h_tN, 0, None)                      # non-negativity

    return h_tN


def fermi_impulse(t, F, t0, k, tau):
    # h(t) = F * exp(-(t - t0)/tau) / (1 + exp((t - t0)/k))
    t = np.asarray(t, dtype=np.float64)
    decay   = np.exp(-np.clip(t - t0, 0, None) / (tau + 1e-8))
    logistic= 1.0 / (1.0 + np.exp((t - t0) / (k + 1e-8)))
    return F * decay * logistic


def fit_fermi_per_pixel(t, h_col, aif_peak_t,
                        k_frac=0.1, tau_frac=0.2, maxfev=3000):
    """
    Fit F, t0, k, tau for one pixel's impulse response h_col (length T).
    aif_peak_t is a reasonable initial guess for arrival scale.
    """
    if not np.any(h_col > 0):
        return np.nan, (np.nan,)*4

    T = len(t)
    duration = t[-1] - t[0]
    F0   = max(h_col[0], np.max(h_col))
    t0_0 = max(0.0, aif_peak_t * 0.7)  # before AIF peak
    k0   = max(1e-3, k_frac * duration)
    tau0 = max(1e-3, tau_frac * duration)
    p0   = (F0, t0_0, k0, tau0)
    bounds = ([0.0, 0.0, 1e-3, 1e-3], [np.inf, t[-1], t[-1], t[-1]])

    try:
        popt, _ = curve_fit(fermi_impulse, t, h_col, p0=p0, bounds=bounds, maxfev=maxfev)
        F, t0, k, tau = popt
        mbf = fermi_impulse(np.array([0.0]), F, t0, k, tau)[0]  # MBF = h_model(0)
        return mbf, popt
    except Exception:
        # Fallback: h(0) as MBF estimate if fit fails
        return float(h_col[0]), (np.nan, ) * 4


def map_mbf(dcm_series_path, tiff_masks_path, is_vis_check=False):
    """
    Returns MBF map of myocard.

    dcm_series_path: path to DICOM images series of perfusion at a single slice (short-axis);

    tiff_masks_path: path to a tiff file that contains 2 binary masks: 1) Left vetricle; 2) Myocard;

    is_vis_check: creates sanity_check.png for one of MYO pixels if set as True (to check that h(t) was found correctly);

    output: numpy.array of MBF map.
    """

    vol_np = read_perf_dcm(dcm_series_path)
    times_np = get_perf_ts(dcm_series_path)

    mask_bpool, mask_myo = read_roi_masks(tiff_masks_path)

    aif_t = get_aif_curve(vol_np, mask_bpool)
    
    vol_myo = np.where(mask_myo, vol_np, 0)
    dt_s = float(times_np[1] - times_np[0])
    T, H, W = vol_myo.shape

    # Ensure mask is boolean and 2-D
    mask_myo = mask_myo.astype(bool)           # (H, W)

    # Build matrix of myocardial curves only where mask is true
    roi_idx = np.argwhere(mask_myo)
    roi_pix_num = roi_idx.shape[0]

    # Extract only masked pixels across time
    T, H, W = vol_myo.shape
    Y_tN = vol_myo[:, mask_myo]                # result: (T, N)
    print(aif_t.shape)                         # (58,)
    print(Y_tN.shape)                          # (58, N)  — good

    # Deconvolution (vectorized)
    alpha = 1e-2
    h_tN = wiener_deconv_all(aif_t, Y_tN, dt_s, alpha=alpha)  # (T, N)

    vol_myo = np.where(mask_myo, vol_np, 0)

    # Initial guess helper: AIF peak time
    aif_peak_t = times_np[int(np.argmax(aif_t))]

    # Fit Fermi per pixel and assemble MBF map
    mbf_map = np.full((H, W), np.nan, dtype=np.float32)
    rho = 1.05  # g/mL

    for col, (y, x) in enumerate(roi_idx):
        h_col = h_tN[:, col]
        mbf_h0_like, popt = fit_fermi_per_pixel(times_np, h_col, aif_peak_t)
        # Convert to mL/g/min
        mbf_map[y, x] = (mbf_h0_like * 60.0) / rho

    if is_vis_check:
        # Sanity check on a random ROI pixel
        idx = np.random.randint(0, roi_pix_num)
        h_one = h_tN[:, idx]
        y_hat = np.convolve(aif_t, h_one)[ : T] * dt_s

        import matplotlib.pyplot as plt
        plt.figure(figsize=(9, 3))
        plt.subplot(1, 3, 1); plt.plot(times_np, aif_t); plt.title("AIF(t)"); plt.xlabel("s")
        plt.subplot(1, 3, 2); plt.plot(times_np, Y_tN[:, idx], label="MYO"); plt.plot(times_np, y_hat, '--', label="recon"); plt.legend(); plt.title("Fit")
        plt.subplot(1, 3, 3); plt.plot(times_np, h_one); plt.title("h(t)"); plt.xlabel("s")
        plt.tight_layout(); 
        plt.savefig("sanity_check.png")

    return mbf_map


if __name__ == "__main__":
    tiff_masks_path = Path("./data/AIF_And_Myo_Masks.tiff")
    dcm_series_path = Path("./data/MotionCorrectedPerfusionSeries/")
    is_vis = True

    mbf_map_np = map_mbf(dcm_series_path, tiff_masks_path, is_vis_check=True)

    if is_vis:
        cv2.imwrite("mbf_map.png", mbf_map_np)

