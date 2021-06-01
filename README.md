# Sentinel-1 SAR Backscatter Analysis Ready Data Preparation in Google Earth Engine

## Introduction
The Sentinel-1 satellites provide temporally dense and high spatial resolution synthetic aperture radar (SAR) imagery. The open data policy and global coverage of Sentinel-1 make it a valuable data source for a wide range of SAR-based applications. In this regard, Google Earth Engine (GEE) is a key platform for large area analysis with preprocessed Sentinel-1 backscatter images being available within few days after acquisition.  In this implementation, we present a framework for preparing Sentinel-1 SAR backscatter Analysis-Ready-Data (ARD) in GEE that implements additional border noise correction, speckle filtering and radiometric terrain normalization. The proposed framework can be used to generate Sentinel-1 ARD suitable for a wide range of land and inland water applications. The ARD preparation framework is implemented in [GEE JavaScript](https://code.earthengine.google.com/?accept_repo=users/adugnagirma/gee_s1_ard) and Python API's.

This framework is intended for researchers and non-experts in microwave remote sensing. It is intended to provide flexibility for a wide variety of large area land and inland water monitoring applications.


## Features
This framework generates a Sentinel-1 SAR ARD by applying three processing modules.
1. Addtional Border noise correction
2. Speckle Filtering 
   - Mono-temporal 
   - Multi-temporal 
3. Radiometric Terrain Normalization

The framework processes single (VV or VH) or dual (VV and VH) polarization data in either ascending, descending or both orbits at the same time. Results can be displayed and exported in the linear or dB scale. 


![flowchart3](https://user-images.githubusercontent.com/48068921/117692979-d840e900-b1bd-11eb-8dd4-a1d552071362.png)

## Usage
The details about parameter setting and their associated methods is described in the main script and accompanying [technical note](https://www.mdpi.com/2072-4292/13/10/1954/htm) published in MDPI Remote sensing. 

To use the framework in GEE code editor, go to the [gee_s1_ard public repo](https://code.earthengine.google.com/?accept_repo=users/adugnagirma/gee_s1_ard) and copy the contents of s1_ard.js to your own repository. The path to the preprocessing functions i.e. ('users/adugnagirma/gee_s1_ard') is a public so you don't need to have the preprocessing functions copied to your repository. 

When using the Python API, the user should adjust the script path and GEE id to their own path and id before processing.

![github_pic2](https://user-images.githubusercontent.com/48068921/117958586-75fdfa80-b31b-11eb-9000-d1eed1ebb675.png)

RGB visualization of a dual polarized (VV and VH) Sentinel-1 SAR backscatter image of central Borneo, Indonesia (Lat: -0.35, Lon: 112.15) (a) as ingested into Google Earth Engine; and (b) after applying additional boarder noise removal, a 9×9 multi-temporal Gamma MAP specklefilter and radiometric terrain normalization with a volume scattering model. Here VV is in red,VH is in green and VV/VH ratio is in blue.

## Dependencies
The JavaScript code runs in the GEE code editor with out installing additional packages. However, the python code requires the installation of 
 [Google Earth Engine](https://github.com/google/earthengine-api) API

## Citation

Mullissa, A.; Vollrath, A.; Odongo-Braun, C.; Slagter, B.; Balling, J.; Gou, Y.; Gorelick, N.; Reiche, J. Sentinel-1 SAR Backscatter Analysis Ready Data Preparation in Google Earth Engine. Remote Sens. 2021, 13, 1954. https://doi.org/10.3390/rs13101954 
