#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar  6 10:55:16 2021

@author: adugnamullissa
"""

import ee
import border_noise_correction as bnc
import speckle_filter as sf
import terrain_flattening as trf

ee.Initialize()

#---------------------------------------------------------------------------//
# Linear to db scale
#---------------------------------------------------------------------------//

def lin_to_db(image):
  """ Convert backscatter from linear to dB. """
  bandNames = image.bandNames().remove('angle')
  db = ee.Image.constant(10).multiply(image.log10()).rename(bandNames)
  return image.addBands(db, None, True)

def lin_to_db2(image):
    """ Convert backscatter from linear to dB by removing the ratio band. """
    if (image.bandNames().contains('VH')):
        db = ee.Image(10).multiply(image.select(['VV', 'VH']).log10()).rename(['VV', 'VH'])
    return image.addBands(db)

def add_ratio_lin(image):
    """ Adding ratio band for visualization """
    ratio = image.addBands(image.select('VV').divide(image.select('VH')).rename('VVVH_ratio'))
    return ratio.set('system:time_start', image.get('system:time_start'))


###########################################
# 1. DO THE JOB
###########################################

def s1_preproc(params):
    """
    Applies preprocessing to a collection of S1 images to return an analysis ready sentinel-1 data. 

    """
    
    APPLY_BORDER_NOISE_CORRECTION = params['APPLY_BORDER_NOISE_CORRECTION']
    APPLY_TERRAIN_FLATTENING = params['APPLY_TERRAIN_FLATTENING']
    APPLY_SPECKLE_FILTERING = params['APPLY_SPECKLE_FILTERING']
    POLARIZATION = params['POLARIZATION']
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
    CLIP_TO_ROI = params['CLIP_TO_ROI']
    ASSET_ID = params['ASSET_ID']

    ###########################################
    # 0. CHECK PARAMETERS
    ###########################################
    
    if APPLY_BORDER_NOISE_CORRECTION is None: APPLY_BORDER_NOISE_CORRECTION = True    
    if APPLY_TERRAIN_FLATTENING is None: APPLY_TERRAIN_FLATTENING = True       
    if APPLY_SPECKLE_FILTERING is None: APPLY_SPECKLE_FILTERING = True         
    if POLARIZATION is None: POLARIZATION = 'VVVH'
    if SPECKLE_FILTER_FRAMEWORK is None: SPECKLE_FILTER_FRAMEWORK = 'MULTI BOXCAR' 
    if SPECKLE_FILTER is None: SPECKLE_FILTER = 'GAMMA MAP' 
    if SPECKLE_FILTER_KERNEL_SIZE is None: SPECKLE_FILTER_KERNEL_SIZE = 7
    if NR_OF_IMAGES is None: NR_OF_IMAGES = 10
    if TERRAIN_FLATTENING_MODEL is None: TERRAIN_FLATTENING_MODEL = 'VOLUME' 
    if TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER is None: TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER = 0 
    if FORMAT is None: FORMAT = 'DB' 
    if ORBIT is None: ORBIT = 'DESCENDING' 
    
    pol_required = ['VV', 'VH', 'VVVH']
    if (POLARIZATION not in pol_required):
        raise ValueError("ERROR!!! Parameter POLARIZATION not correctly defined")

    model_required = ['DIRECT', 'VOLUME']
    if (TERRAIN_FLATTENING_MODEL not in model_required):
        raise ValueError("ERROR!!! Parameter TERRAIN_FLATTENING_MODEL not correctly defined")

    format_required = ['LINEAR', 'DB']
    if (FORMAT not in format_required):
        raise ValueError("ERROR!!! FORMAT not correctly defined")
        
    frame_needed = ['MONO', 'MULTI']
    if (SPECKLE_FILTER_FRAMEWORK not in frame_needed.keys):
        raise ValueError("ERROR!!! SPECKLE_FILTER_FRAMEWORK not correctly defined")

    format_sfilter = ['MONO BOXCAR', 'MONO LEE', 'MONO GAMMA MAP' \
              ,'MONO REFINED LEE', 'MONO LEE SIGMA', 'MULTI BOXCAR', 'MULTI LEE' \
              ,'MULTI GAMMA MAP', 'MULTI REFINED LEE' \
              ,'MULTI LEE SIGMA']
    if (SPECKLE_FILTER not in format_sfilter):
        raise ValueError("ERROR!!! SPECKLE_FILTER not correctly defined")

    if (TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER < 0):
        raise ValueError("ERROR!!! TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER not correctly defined")

    if (SPECKLE_FILTER_KERNEL_SIZE <= 0):
        raise ValueError("ERROR!!! SPECKLE_FILTER_KERNEL_SIZE not correctly defined")
    
    
    ###########################################
    # 1. IMPORT COLLECTION
    ###########################################
    
    s1 = ee.ImageCollection('COPERNICUS/S1_GRD_FLOAT') \
                .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH')) \
                .filter(ee.Filter.eq('orbitProperties_pass', ORBIT)) \
                .filter(ee.Filter.eq('instrumentMode', 'IW')) \
                .filter(ee.Filter.eq('resolution_meters', 10)) \
                .filterDate(START_DATE, STOP_DATE) \
                .filterBounds(ROI)
    
    #########################
    # 2. SELECT POLARIZATION
    #########################
    
    if (POLARIZATION == 'VV'):
      s1 = s1.select(['VV','angle'])
    elif (POLARIZATION == 'VH'):
      s1 = s1.select(['VH','angle'])
    elif (POLARIZATION == 'VVVH'):
      s1 = s1.select(['VV','VH','angle'])  
      
    ########################
    # 3. CLIP TO ROI
    ####################### 
    
    if (CLIP_TO_ROI == True):
      s1 = s1.map(lambda image: image.clip(ROI))
    
    ###########################################
    # 4. ADDITIONAL BORDER NOISE CORRECTION
    ###########################################
    
    if (APPLY_BORDER_NOISE_CORRECTION == True):
      s1_1 = s1.map(bnc.f_mask_edges)
    else : 
      s1_1 = s1
          
      
    ########################
    # 5. SPECKLE FILTERING
    ####################### 
    
    if (APPLY_SPECKLE_FILTERING) :
        if (SPECKLE_FILTER_FRAMEWORK == 'MONO') :
            s1_1 = ee.ImageCollection(sf.MonoTemporal_Filter(s1_1, SPECKLE_FILTER_KERNEL_SIZE, SPECKLE_FILTER ))
        else :
            s1_1 = ee.ImageCollection(sf.MultiTemporal_Filter(s1_1, SPECKLE_FILTER_KERNEL_SIZE, SPECKLE_FILTER,NR_OF_IMAGES ));  

              
    ########################
    # 6. TERRAIN CORRECTION
    ####################### 
    
    if (APPLY_TERRAIN_FLATTENING == True):
        s1_1 = trf.slope_correction(s1_1 \
                              , TERRAIN_FLATTENING_MODEL \
                              ,DEM \
                              ,POLARIZATION \
                              ,TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER)
              
    ########################
    # 7. FORMAT CONVERSION
    ####################### 
    
    if (FORMAT == 'DB'):
        s1_1 = s1_1.map(lin_to_db)
        
    ########################
    # 8. EXPORT ASSET
    #######################  
    

    if (SAVE_ASSET == True): 
            
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
        