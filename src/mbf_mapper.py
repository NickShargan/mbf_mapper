import argparse
from pathlib import Path

import numpy as np
from numpy.fft import rfft, irfft
from scipy.optimize import curve_fit

from utils import read_perf_dcm, get_perf_ts, read_roi_masks
from vis import plot_mbf_map, plot_sanity_check


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
    """Fermi function defined in Jerosch-Herold, 1998"""
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
        y_recon = np.convolve(aif_t, h_one)[ : T] * dt_s

        plot_sanity_check(times_np, aif_t, Y_tN, idx, y_recon, h_one)

    return mbf_map


def main():
    parser = argparse.ArgumentParser(
        description="Compute myocardial blood flow (MBF) map from DICOM perfusion series and TIFF masks."
    )

    parser.add_argument(
        "--dcm_series",
        type=Path,
        required=True,
        help="Path to DICOM perfusion series directory."
    )
    parser.add_argument(
        "--tiff_masks",
        type=Path,
        required=True,
        help="Path to TIFF file containing left ventricle and myocardium masks."
    )
    parser.add_argument(
        "--vis_check",
        action="store_true",
        help="If set, saves sanity_check.png for visualization."
    )
    # parser.add_argument(
    #     "--alpha",
    #     type=float,
    #     default=1e-2,
    #     help="Wiener deconvolution regularization parameter (default: 1e-2)."
    # )
    parser.add_argument(
        "--out_dir",
        type=Path,
        default=Path("./"),
        help="Output directory for results (default: current directory)."
    )

    args = parser.parse_args()

    mbf_map_np = map_mbf(
        args.dcm_series,
        args.tiff_masks,
        is_vis_check=args.vis_check
    )

    if args.vis_check:
        out_map_path = args.out_dir / "mbf_map.png"
        plot_mbf_map(mbf_map_np, out_map_path)
        print(f"Saved MBF map to {out_map_path}")


if __name__ == "__main__":
    main()
