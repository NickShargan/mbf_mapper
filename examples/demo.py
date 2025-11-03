from pathlib import Path

from mbf_mapper import map_mbf, plot_mbf_map


def main():

    dcm_series_fpath = Path("./data/MotionCorrectedPerfusionSeries/")
    tiff_masks_fpath = Path("./data/AIF_And_Myo_Masks.tiff")

    mbf_map_np = map_mbf(dcm_series_fpath, tiff_masks_fpath, is_vis_check=True)

    out_map_path = "mbf_map.png"
    plot_mbf_map(mbf_map_np, out_map_path)
    print(f"Saved MBF map to {out_map_path}")


if __name__ == "__main__":
    main()
