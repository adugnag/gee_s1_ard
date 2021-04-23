#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Version: v1.2
Date: 2021-04-01
Authors: Mullissa A., Vollrath A., Gorelick N.,  Reiche J., Slagter B., Balling J. , Gou Y., Braun, C.
Description: A wrapper function to derive the Sentinel-1 ARD
"""

import ee
import border_noise_correction as bnc
import speckle_filter as sf
import terrain_flattening as trf
import helper

ee.Initialize()


###########################################
# DO THE JOB
###########################################

def s1_preproc(params):
    """
    Applies preprocessing to a collection of S1 images to return an analysis ready sentinel-1 data.

    Parameters
    ----------
    params : List
        These parameters determine the data selection and image processing parameters.

    Raises
    ------
    ValueError
        

    Returns
    -------
    ee.ImageCollection
        A processed Sentinel-1 image collection

    """
    
    APPLY_BORDER_NOISE_CORRECTION = params['APPLY_BORDER_NOISE_CORRECTION']
    APPLY_TERRAIN_FLATTENING = params['APPLY_TERRAIN_FLATTENING']
    APPLY_SPECKLE_FILTERING = params['APPLY_SPECKLE_FILTERING']
    POLARIZATION = params['POLARIZATION']
    ORBIT = params['ORBIT']
    SPECKLE_FILTER_FRAMEWORK = params['SPECKLE_FILTER_FRAMEWORK']
    SPECKLE_FILTER = params['SPECKLE_FILTER']
    SPECKLE_FILTER_KERNEL_SIZE = params['SPECKLE_FILTER_KERNEL_SIZE']
    NR_OF_IMAGES = params['NR_OF_IMAGES']
    TERRAIN_FLATTENING_MODEL = params['TERRAIN_FLATTENING_MODEL']
    DEM = params['DEM']
    TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER = params['TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER']
    FORMAT = params['FORMAT']
    START_DATE = params['START_DATE']
    STOP_DATE = params['STOP_DATE']
    ORBIT = params['ORBIT']
    ROI = params['ROI']
    SAVE_ASSET = params['SAVE_ASSET']
    ASSET_ID = params['ASSET_ID']

    ###########################################
    # 0. CHECK PARAMETERS
    ###########################################

    if APPLY_BORDER_NOISE_CORRECTION is None:
        APPLY_BORDER_NOISE_CORRECTION = True
    if APPLY_TERRAIN_FLATTENING is None:
        APPLY_TERRAIN_FLATTENING = True
    if APPLY_SPECKLE_FILTERING is None:
        APPLY_SPECKLE_FILTERING = True
    if POLARIZATION is None:
        POLARIZATION = 'VVVH'
    if ORBIT is None:
        ORBIT = 'BOTH'
    if SPECKLE_FILTER_FRAMEWORK is None:
        SPECKLE_FILTER_FRAMEWORK = 'MULTI BOXCAR'
    if SPECKLE_FILTER is None:
        SPECKLE_FILTER = 'GAMMA MAP'
    if SPECKLE_FILTER_KERNEL_SIZE is None:
        SPECKLE_FILTER_KERNEL_SIZE = 7
    if NR_OF_IMAGES is None:
        NR_OF_IMAGES = 10
    if TERRAIN_FLATTENING_MODEL is None:
        TERRAIN_FLATTENING_MODEL = 'VOLUME'
    if TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER is None:
        TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER = 0
    if FORMAT is None:
        FORMAT = 'DB'
    if ORBIT is None:
        ORBIT = 'DESCENDING'

    pol_required = ['VV', 'VH', 'VVVH']
    if (POLARIZATION not in pol_required):
        raise ValueError("ERROR!!! Parameter POLARIZATION not correctly defined")

    orbit_required = ['ASCENDING', 'DESCENDING', 'BOTH']
    if (ORBIT not in orbit_required):
        raise ValueError("ERROR!!! Parameter ORBIT not correctly defined")

    model_required = ['DIRECT', 'VOLUME']
    if (TERRAIN_FLATTENING_MODEL not in model_required):
        raise ValueError("ERROR!!! Parameter TERRAIN_FLATTENING_MODEL not correctly defined")

    format_required = ['LINEAR', 'DB']
    if (FORMAT not in format_required):
        raise ValueError("ERROR!!! FORMAT not correctly defined")

    frame_needed = ['MONO', 'MULTI']
    if (SPECKLE_FILTER_FRAMEWORK not in frame_needed):
        raise ValueError("ERROR!!! SPECKLE_FILTER_FRAMEWORK not correctly defined")

    format_sfilter = ['MONO BOXCAR', 'MONO LEE', 'MONO GAMMA MAP'
              ,'MONO REFINED LEE', 'MONO LEE SIGMA', 'MULTI BOXCAR', 'MULTI LEE'
              ,'MULTI GAMMA MAP', 'MULTI REFINED LEE'
              ,'MULTI LEE SIGMA']
    if (SPECKLE_FILTER not in format_sfilter):
        raise ValueError("ERROR!!! SPECKLE_FILTER not correctly defined")

    if (TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER < 0):
        raise ValueError("ERROR!!! TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER not correctly defined")

    if (SPECKLE_FILTER_KERNEL_SIZE <= 0):
        raise ValueError("ERROR!!! SPECKLE_FILTER_KERNEL_SIZE not correctly defined")

    ###########################################
    # 1. DATA SELECTION
    ###########################################

    # select S-1 image collection
    s1 = ee.ImageCollection('COPERNICUS/S1_GRD_FLOAT') \
        .filter(ee.Filter.eq('instrumentMode', 'IW')) \
        .filter(ee.Filter.eq('resolution_meters', 10)) \
        .filterDate(START_DATE, STOP_DATE) \
        .filterBounds(ROI)

    # select orbit
    if (ORBIT != 'BOTH'):
        s1 = s1.filter(ee.Filter.eq('orbitProperties_pass', ORBIT))

    # select polarization
    if (POLARIZATION == 'VV'):
        s1 = s1.select(['VV', 'angle'])
    elif (POLARIZATION == 'VH'):
        s1 = s1.select(['VH', 'angle'])
    elif (POLARIZATION == 'VVVH'):
        s1 = s1.select(['VV', 'VH', 'angle'])

    ###########################################
    # 2. ADDITIONAL BORDER NOISE CORRECTION
    ###########################################

    if (APPLY_BORDER_NOISE_CORRECTION):
        s1_1 = s1.map(bnc.f_mask_edges)
    else:
        s1_1 = s1
    ########################
    # 3. SPECKLE FILTERING
    #######################

    if (APPLY_SPECKLE_FILTERING):
        if (SPECKLE_FILTER_FRAMEWORK == 'MONO'):
            s1_1 = ee.ImageCollection(sf.MonoTemporal_Filter(s1_1, SPECKLE_FILTER_KERNEL_SIZE, SPECKLE_FILTER))
        else:
            s1_1 = ee.ImageCollection(sf.MultiTemporal_Filter(s1_1, SPECKLE_FILTER_KERNEL_SIZE, SPECKLE_FILTER, NR_OF_IMAGES))

    ########################
    # 4. TERRAIN CORRECTION
    #######################

    if (APPLY_TERRAIN_FLATTENING):
        s1_1 = trf.slope_correction(s1_1 
                                    ,TERRAIN_FLATTENING_MODEL
                                        ,DEM
                                                ,TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER)

    ########################
    # 5. OUTPUT
    #######################

    if (FORMAT == 'DB'):
        s1_1 = s1_1.map(helper.lin_to_db)

    if (SAVE_ASSET):
        size = s1_1.size()
        imlist = s1_1.toList(size)

        for idx in range(0, size):
            img = imlist.get(idx)
            img = ee.Image(img)
            name = img.id()
            description = name
            assetId = ASSET_ID+"/"+name

            task = ee.batch.Export.image.toAsset(image=img,
                                                 assetId=assetId,
                                                 description=description,
                                                 region=img.geometry(),
                                                 scale=10,
                                                 maxPixels=1e13)
            task.start()
            print('Exporting {} to {}'.format(name, assetId))
    return s1_1
