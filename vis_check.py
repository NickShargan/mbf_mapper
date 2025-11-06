from pathlib import Path
import pandas as pd
import numpy as np
import SimpleITK as sitk

# from .src.vis import plot_sanity_check

import plotly.graph_objects as go
from plotly.subplots import make_subplots


def plot_sanity_check(times_np, aif_t, myo_pix_t, y_recon, h_one):
    """Creates sanity_check.png for a random MYO pixel to ensure that h(t) was found correctly."""
    
    fig = make_subplots(
        rows=1, cols=3,
        subplot_titles=["AIF(t)", "Fit", "h(t)"],
        horizontal_spacing=0.08,
    )

    # 1) AIF(t)
    fig.add_trace(
        go.Scatter(x=times_np, y=aif_t, mode="lines", name="AIF"),
        row=1, col=1
    )

    # 2) Fit: measured vs recon
    fig.add_trace(
        go.Scatter(x=times_np, y=myo_pix_t, mode="lines", name="MYO"),
        row=1, col=2
    )
    fig.add_trace(
        go.Scatter(x=times_np, y=y_recon, mode="lines", name="recon", line=dict(dash="dash")),
        row=1, col=2
    )

    # 3) h(t)
    fig.add_trace(
        go.Scatter(x=times_np, y=h_one, mode="lines", name="h(t)"),
        row=1, col=3
    )

    # axes labels + layout
    fig.update_xaxes(title_text="s", row=1, col=1)
    fig.update_xaxes(title_text="s", row=1, col=2)
    fig.update_xaxes(title_text="s", row=1, col=3)

    fig.update_layout(
        height=350,
        width=1000,
        showlegend=True,
        title_text="Sanity check on a random ROI pixel",
    )

    fig.write_image("sanity_check.png", scale=2)
    print("saved to sanity_check.png")


def plot_mbf_map(mbf_map, out_map_path=Path("mbf_map.png")):
    """Creates mbf_map.png with defined MBF map."""
    fig = go.Figure(
        data=go.Heatmap(
            z=np.flipud(mbf_map),
            colorscale="Hot",
            colorbar=dict(title="mL/g/min"),
        )
    )

    fig.update_layout(
        title="MBF (mL/g/min)",
        height=500,
        width=500,
        # xaxis=dict(showticklabels=False),
        # yaxis=dict(showticklabels=False),
    )

    fig.write_image(out_map_path, scale=2)


def main():
    aif_cpp_df = pd.read_csv("./build/aif.csv")
    myo_pix_cpp_df = pd.read_csv("./build/myo_pix.csv")
    h_t_cpp_df = pd.read_csv("./build/h_t.csv")

    times = aif_cpp_df["Time_s"].values
    aif_t = aif_cpp_df["AIF"].values
    myo_pix_t = myo_pix_cpp_df["MYO_pix"].values
    h_t = h_t_cpp_df["h_t_svd"].values

    myo_pix_conv_t = np.convolve(aif_t, h_t)

    plot_sanity_check(times, aif_t, myo_pix_t, myo_pix_conv_t, h_t)

    img = sitk.ReadImage("./cpp/mbf_map.nii")
    mbf_map = sitk.GetArrayFromImage(img)
    plot_mbf_map(mbf_map)


if __name__ == "__main__":
    main()
