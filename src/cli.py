# src/cli.py
import argparse
import sys
from pathlib import Path

from mbf_mapper import map_mbf
from vis import plot_mbf_map


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

    return 0


if __name__ == "__main__":
    sys.exit(main())