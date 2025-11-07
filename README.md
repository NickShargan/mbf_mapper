# MBF mapper

This package can be used to estimate myocardial perfusion and visualise corresponding map:


<!-- ### Input

Perfusion DICOM series, blood pool (e.g. left ventricle) & Myocardium segmentation masks:

![image](./docs/perfusion.gif) ![image](./docs/ventr_mask.png) ![image](./docs/myo_mask.png)

Blood pool (e.g. left ventricle) & Myocardium segmentation masks: -->

## Compile

Build and run:

```
mkdir build
cd build
cmake -DITK_DIR=~/ITK/build ..          
make 
```

## Usage
Ensure that MotionCorrectedPerfusionSeries/ and AIF_And_Myo_Masks.tiff are placed into ./data/ directory. Then, run:

```
./mbf_cli ../data/MotionCorrectedPerfusionSeries/ ../data/AIF_And_Myo_Masks.tiff ../mbf_map.nii
```

### Sanity check for one of pixels

```
# install dependencies if needed
python vis_check.py
```

h(t) - was defined by SVD deconvolution; impulse response function

MYO - concentration of contrastive substance at a specific pixel/voxel

recon - reconstruction of MYO obtained by convoluting AIF and h(t)

![image](./docs/sanity_check_single_pix.png)

### MBF calculation (in ml/min/gr)
MBF for this particular pixel was calculated based on h(0) value by the following logic:

```
# tissue density:
const double rho = 1.05

# dt - time difference between neighbouring images/slices (0th and 1st)

double mbf = (h0 / dt) * 60.0 / rho;
```

### MBF map

The defined MBF map is written to mbf_map.nii . It can be displayed in Slicer3D, Python(SimpleITK & Plotly) or by other visualisational tools.

![image](./docs/mbf_map.png)
