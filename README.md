# QLFM
matlab code for **"Computational optical sectioning with an incoherent multiscale scattering model for light-field microscopy"**
Version: 1.0 Copyright: 2021, YI ZHANG, ZHI LU, JIAMINWU, XING LIN, DONG JIANG, YEYI CAI, JIACHEN XIE, YULING WANG, TIANYI ZHU, XIANGYANG JI, QIONGHAI DAI

## Contents
- [Overview](#Overview)
- [Directory structure](#Directory-structure)
- [How to use](#How-to-use)
- [Citation](#Citation)

## Overview
This repository provides the implementation of the QLFM reconstruction with different models, including QLFM without scattering using ideal PSFs, QLFM without scattering using calibrated PSFs, QLFM with alternating scattering directions.

Additional example data and experimental PSFs can be downloaded with the following link:
https://drive.google.com/drive/folders/1bTH5eibetRsKuxwTPZk9sAXzOvdA2wgf?usp=sharing

## Directory structure

```
QLFM
|---recon_main
|---|---backwardProjectACC_fft2.m
|---|---deconvGPU_QLFMScatter_2D.m
|---|---deconvGPU_QLFMScatter_L2R.m
|---|---deconvGPU_QLFMScatter_R2L.m
|---|---deconvGPU_QLFMwithoutScatter.m
|---|---forwardProjectACC_fft2.m
|---|---imwrite3dTIFF.m
|---|---main.m
|---|---map.mat
|---|---SampleV.m
|---|---Recon_QLFMwithoutScatter_cali
|---|---|---c1_frame0.tif
|---|---Recon_QLFMwithoutScatter_ideal
|---|---|---c1_frame0.tif
|---|---Recon_QLFMwithScatter_2directions
|---|---|---c1_frame0.tif
|---|---Recon_QLFMwithScatter_Left2Right
|---|---|---c1_frame0.tif
|---|---Recon_QLFMwithScatter_Right2Left
|---|---|---c1_frame0.tif
|---psf
|---|---index.mat
|---|---psf_cali2.mat
|---|---psf_ideal2.mat
|---|---results
|---dataset
|---|---c1-1.tif
|---README.md
```
## How to use
### Environment
*   MATLAB 2020b (64bit)
*   Windows 10 64bit 
*   Intel i9-9820X CPU
*   NVIDIA TITAN RTX GPU
*   128G RAM
## Citation
If you use the code, please cite the xxx
