# Sentinel-1 SAR Backscatter Analysis-Ready-Data Preparation in Google Earth Engine

## Introduction
The Sentinel-1 satellites provide temporally dense and high spatial resolution synthetic aperture radar (SAR) imagery. The open data policy and global coverage of Sentinel-1 make it a valuable data source for a wide range of SAR-based applications. In this regard, Google Earth Engine (GEE) is a key platform for large area analysis with preprocessed Sentinel-1 backscatter images being available within few days after acquisition.  In this implementation, we present a framework for preparing Sentinel-1 SAR backscatter Analysis-Ready-Data in Google Earth Engine that implements additional border noise correction, speckle filtering and radiometric terrain normalization. The proposed framework can be used to generate Sentinel-1 Analysis-Ready-Data suitable for a wide range of land and inland water applications. The ARD preparation framework is implemented in Google Earth Engine JavaScript and Python API's.

This framework is intended for students and researchers, who are not experts in microwave remote sensing. It may also be useful to microwaave remote sensing experts who are experienced in GEE to use this framework as afoundation for advanced implementations.

## Features
This framework consists of three processing modules.
1. Addtional Border noise correction
2. Speckle Filtering
    i. Mono-temporal speckle filters
       -Boxcar
       -Lee
       -Gamma MAP
       -Refined Lee
       -Improved Lee sigma
    ii. Multi-temporal speckle filters
       -Quegan and Yu
3. Radiometric Terrain Normalization

## Usage
The Following parameters should be filled as discussed below.
START_DATE: The earliest date to include images for (inclusive).<br/>
END_DATE: The latest date to include images for (exclusive).<br/>
POLARIZATION: The Sentinel-1 image polarization to select for processing.<br/>
            'VV' - selects the VV polarization.<br/>
            'VH' - selects the VH polarization.<br/>
            "VVVH' - selects both the VV and VH polarization for processing.<br/>
ORBIT:  The orbits to include. (string: BOTH, ASCENDING or DESCENDING)<br/>
GEOMETRY: The region to include imagery within.<br/>
            The user can interactively draw a bounding box within the map window or define the edge coordinates.<br/>
APPLY_BORDER_NOISE_CORRECTION: (Optional) true or false options to apply additional Border noise correction:<br/>
APPLY_SPECKLE_FILTERING: (Optional) true or false options to apply speckle filter<br/>
SPECKLE_FILTER: Type of speckle filtering to apply (String). If the APPLY_SPECKLE_FILTERING parameter is true then the selected speckle filter type will be used.<br/>
            'BOXCAR' - Applies a boxcar filter on each individual image in the collection<br/>
            'LEE' - Applies a Lee filter on each individual image in the collection based on <br/>
            'GAMMA MAP' - Applies a Gamma maximum a-posterior speckle filter on each individual image in the collection <br/>
            'REFINED LEE' - Applies the Refined Lee speckle filter on each individual image in the collection<br/>
            'LEE SIGMA' - Applies the improved Lee sigma speckle filter on each individual image in the collection<br/>
SPECKLE_FILTER_FRAMEWORK: is the framework where filtering is applied (String). It can be 'MONO' or 'MULTI'. In the MONO case
                          the filtering is applied to each image in the collection individually. Whereas, in the MULTI case,
                          the Multitemporal Speckle filter is applied based on  [6] with any of the above mentioned speckle filters.
SPECKLE_FILTER_KERNEL_SIZE: is the size of the filter spatial window applied in speckle filtering. It must be a positive odd integer.<br/>
NR_OF_IMAGES: is the number of images to use in the multi-temporal speckle filter framework.<br/>
TERRAIN_FLATTENING : (Optional) true or false option to apply Terrain correction.
TERRAIN_FLATTENING_MODEL : model to use for radiometric terrain normalization (DIRECT, or VOLUME)<br/>
DEM : digital elevation model (DEM) to use (as EE asset)<br/>
TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER : additional buffer parameter for passive layover/shadow mask in meters<br/>
FORMAT : the output format for the processed collection. this can be 'LINEAR' or 'DB'.<br/>
CLIP_TO_ROI: (Optional) Clip the processed image to the region of interest.<br/>
SAVE_ASSETS : (Optional) Exports the processed collection to an asset.<br/>
ASSET_ID : (Optional) The user id path to save the assets.<br/>
        
The processing returns an ee.ImageCollection with an analysis ready Sentinel 1 imagery with the specified polarization images and angle band.

## Dependencies
The JavaScript code runs in the GEE code editor with out installing additional packages. However, the python code requires <br/>
 [Google Earth Engine](https://github.com/google/earthengine-api)

## References
If you find this implementation useful please cite our work as

Mullissa, A., Vollrath, A., Odongo-Braun, C., Slagter, B., Balling, J., Gou, Y., Gorelick, N., Reiche, J. (2021). Sentinel-1 SAR Backscatter Analysis-Ready-Data Preparation inGoogle Earth Engine. Remote Sensing, xxx, 1-7.
