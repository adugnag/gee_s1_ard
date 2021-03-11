#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar  6 10:55:17 2021

@author: adugnamullissa

    Applies preprocessing to a collection of S1 images to return an analysis ready sentinel-1 data. 

    params:
        ROI: The region to include imagery within. (ee.Geometry) The user has to interactively draw a bounding box within the map window.
        ORBITS:  The orbits to include. (string: ASCENDING or DESCENDING)
        START_DATE: The earliest date to include images for (inclusive).
        END_DATE: The latest date to include images for (exclusive).
        POLARIZATION: The Sentinel-1 image polarization to select for processing.
            'VV' - selects the VV polarization.
            'VH' - selects the VH polarization.
            "VVVH' - selects both the VV and VH polarization for processing.
        APPLY_SPECKLE_FILTERING: true or false options to apply speckle filter
        SPECKLE_FILTER: (Optional)   Type of speckle filtering to apply. Can be monotemporal or multitemporal.
            This is done by passing the argument 'MONO' or 'MULTI'.
            'MONO BOXCAR' - implements a boxcar filter on each individual image in the collection
            'MONO LEE' - implements a Lee filter on each individual image in the collection based on (J.S. Lee et al. 1980)
            'MONO GAMMA MAP' - implements a Gamma maximum a-posterior speckle filter on each individual image in the collection based on (Lopez et al. 1990
            'MONO REFINED LEE' - implements the Refined Lee speckle filter on each individual image in the collection
                                  based on (J.S.Lee et al. 1999)
            'MONO LEE SIGMA' - implements the improved Lee sigma speckle filter on each individual image in the collection
                                  based on (J.S.Lee et al. 2009)
            'MULTI' - Multitemporal Speckle filter based on the approach by Quegan and Yu 2001 with any of the above mentioned speckle filters.
        SPECKLE_FILTER_KERNEL_SIZE: is the size of the filter spatial window applied in speckle filtering. must be a positive integer.
        APPLY_BORDER_NOISE_CORRECTION: (Optional) true or false options to apply additional Border noise correction:
        TERRAIN_FLATTENING - (Optional) true or false option to apply Terrain correction based on Vollrath et al. 2020 and Hoekman and Reiche 2015 is applied. 
        TERRAIN_FLATTENING_MODEL - which model to use for radiometric terrain flattening (DIRECT, or VOLUME)
        DEM - which digital elevation model (DEM) to use (as EE asset)
        TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER - buffer parameter for layover/shadow mask in meters
        FORMAT - the output format for the processed collection. this can be 'LINEAR' or 'DB'.
        CLIP_TO_ROI - (Optional) Resizes the processed images to the region of interest.
        SAVE_ASSETS - (Optional) Exports the processed collection to an asset.

    Returns:
        (s1_processed) where
            s1_processed (ee.ImageCollection): preprocessed image collection

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
            'CLIP_TO_ROI': False,
            'SAVE_ASSET': False, 
            'ASSET_ID': "users/adugnagirma"
            }
#processed s1 collection
s1_processed = wp.s1_preproc(params)
