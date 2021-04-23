#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Version: v1.2
Date: 2021-03-12
Authors: Mullissa A., Vollrath A.,  Reiche J., Slagter B., Balling J. , Gou Y., Braun, C.
Description: This script creates an analysis ready S1 image collection.

    Params:
        geometry: The region to include imagery within. (ee.Geometry); The user has to interactively draw a bounding box within the map window.
        ORBIT:  The orbits to include. (string:BOTH, ASCENDING or DESCENDING)
        START_DATE: The earliest date to include images for (inclusive).
        END_DATE: The latest date to include images for (exclusive).
        POLARIZATION: The Sentinel-1 image polarization to select for processing.
            'VV' - selects the VV polarization.
            'VH' - selects the VH polarization.
            "VVVH' - selects both the VV and VH polarization for processing.
        APPLY_SPECKLE_FILTERING: true or false options to apply speckle filter
        SPECKLE_FILTER: (Optional)   Type of speckle filtering to apply. 
            'BOXCAR' - implements a boxcar filter on each individual image in the collection
            'LEE' - implements a Lee filter on each individual image in the collection based on (J.S. Lee et al. 1980)
            'GAMMA MAP' - implements a Gamma maximum a-posterior speckle filter on each individual image in the collection based on (Lopez et al. 1990
            'REFINED LEE' - implements the Refined Lee speckle filter on each individual image in the collection
                                  based on (J.S.Lee et al. 1999)
            'LEE SIGMA' - implements the improved Lee sigma speckle filter on each individual image in the collection
                                  based on (J.S.Lee et al. 2009)
        SPECKLE_FILTER_FRAMEWORK: is the framework where filtering is applied. it can be 'MONO' or 'MULTI'. In the MONO case
                                  the filtering is applied to each image in the collection individually. Whereas, the 
                                  Multitemporal Speckle filter based on the approach by Quegan and Yu 2001 with any of the above mentioned speckle filters.
        SPECKLE_FILTER_KERNEL_SIZE: is the size of the filter spatial window applied in speckle filtering. must be a positive integer.
        NR_OF_IMAGES: is the number of images to use in the multi-temporal filter framework.
        APPLY_BORDER_NOISE_CORRECTION: (Optional) true or false options to apply additional Border noise correction:
        TERRAIN_FLATTENING : (Optional) true or false option to apply Terrain correction based on Vollrath et al. 2020 and Hoekman and Reiche 2015 is applied. 
        TERRAIN_FLATTENING_MODEL : which model to use for radiometric terrain flattening (DIRECT, or VOLUME)
        DEM : which digital elevation model (DEM) to use (as EE asset)
        TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER : buffer parameter for layover/shadow mask in meters
        FORMAT : the output format for the processed collection. this can be 'LINEAR' or 'DB'.
        SAVE_ASSETS : (Optional) Exports the processed collection to an asset.
    Returns:
        An ee.ImageCollection with an analysis ready Sentinel 1 imagery with the specified polarization images and angle band.

    """

import wrapper as wp
import ee

ee.Initialize()

#/***************************/ 
#// MAIN
#/***************************/ 
#Parameters
params = {  'START_DATE': '2018-01-01', 
            'STOP_DATE': '2018-01-15',        
            'ORBIT': 'DESCENDING',
            'POLARIZATION': 'VVVH',
            'ORBIT' : 'BOTH',
            'ROI': ee.Geometry.Rectangle([-47.1634, -3.00071, -45.92746, -5.43836]),
            'DEM': ee.Image('USGS/SRTMGL1_003'),
            'APPLY_BORDER_NOISE_CORRECTION': False,
            'APPLY_TERRAIN_FLATTENING': True,
            'TERRAIN_FLATTENING_MODEL': 'VOLUME',
            'TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER':0,
            'APPLY_SPECKLE_FILTERING': True,
            'SPECKLE_FILTER_FRAMEWORK':'MULTI',
            'SPECKLE_FILTER': 'BOXCAR',
            'SPECKLE_FILTER_KERNEL_SIZE': 7,
            'NR_OF_IMAGES':10,
            'FORMAT': 'DB',
            'SAVE_ASSET': False, 
            'ASSET_ID': "users/adugnagirma"
            }
#processed s1 collection
s1_processed = wp.s1_preproc(params)
