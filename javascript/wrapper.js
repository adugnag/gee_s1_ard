// File: S1_ARD_TEST/helper.js
// Version: v1.2
// Date: 2021-02-11
// Authors: Mullissa A., Vollrath A.,  Reiche J., Slagter B., Balling J. , Gou Y., Braun, C.

//****************************
// ALL PREPROCESSING
//****************************

var speckle_filter = require('users/adugnagirma/s1_ard_master:speckle_filter');
var terrain_flattening = require('users/adugnagirma/s1_ard_master:terrain_flattening');
var border_noise_correction = require('users/adugnagirma/s1_ard_master:border_noise_correction');

exports.s1_preproc = function(s1, APPLY_BORDER_NOISE_CORRECTION, 
                                  APPLY_TERRAIN_FLATTENING, 
                                  APPLY_SPECKLE_FILTERING,
                                  SPECKLE_FILTER_FRAMEWORK,
                                  POLARIZATION,
                                  SPECKLE_FILTER,
                                  SPECKLE_FILTER_KERNEL_SIZE,
                                  NR_OF_IMAGES,
                                  TERRAIN_FLATTENING_MODEL, 
                                  DEM,
                                  TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER, 
                                  FORMAT) {
  
  /************************  
   * 0. CHECK PARAMETERS  
  ************************/
  if (SPECKLE_FILTER === undefined) SPECKLE_FILTER = "GAMMA MAP";
  if (SPECKLE_FILTER_KERNEL_SIZE===undefined) SPECKLE_FILTER_KERNEL_SIZE = 7;
  if (TERRAIN_FLATTENING_MODEL === undefined) TERRAIN_FLATTENING_MODEL = 'VOLUME';
  if (TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER === undefined) TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER = 0;
  if (FORMAT === undefined) FORMAT = 'DB';
  if (POLARIZATION === undefined) POLARIZATION = 'VVVH';
  if (APPLY_BORDER_NOISE_CORRECTION === undefined) APPLY_BORDER_NOISE_CORRECTION = true;
  if (APPLY_TERRAIN_FLATTENING===undefined) APPLY_TERRAIN_FLATTTENING = true;
  if (APPLY_SPECKLE_FILTERING===undefined) APPLY_SPECKLE_FILTERING = true; 
  if (SPECKLE_FILTER_FRAMEWORK===undefined) SPECKLE_FILTER_FRAMEWORK = 'MULTI';
  
  function notContains(list, value) {return list.indexOf(value) == -1;}
  
  var pol_required = ['VV', 'VH', 'VVVH']
  if (notContains(pol_required, POLARIZATION)) {
       throw new Error("Parameter POLARIZATION not correctly defined")
  } 
  
  var model_required = ['DIRECT', 'VOLUME']
  if (notContains(model_required, TERRAIN_FLATTENING_MODEL)) {
       throw new Error("Parameter TERRAIN_FLATTENING_MODEL not correctly defined")
  } 
  
  var format_required = ['LINEAR', 'DB']
  if (notContains(format_required, FORMAT)) {
       throw new Error("Parameter FORMAT not correctly defined")
  } 
  
  var frame_required = ['MONO', 'MULTI']
  if (notContains(frame_required, SPECKLE_FILTER_FRAMEWORK)) {
       throw new Error("Parameter SPECKLE_FILTER_FRAMEWORK not correctly defined")
  } 
  
  var sfilter_required = ['BOXCAR', 'LEE', 'GAMMA MAP'
              , 'REFINED LEE', 'LEE SIGMA', 'BOXCAR', 'LEE'
              , 'GAMMA MAP', 'REFINED LEE'
              , 'LEE SIGMA']
  if (notContains(sfilter_required, SPECKLE_FILTER)) {
       throw new Error("Parameter SPECKLE_FILTER not correctly defined")
  } 

  if (TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER < 0) {
  throw new Error("The TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER should be GREATER or EQUAL to 0")}

  if (SPECKLE_FILTER_KERNEL_SIZE <= 0) {
  throw new Error("The SPECKLE_FILTER_KERNEL_SIZE should be a positive integer")}
  
  /************************  
   * 1. SELECT POLARIZATION  
  ************************/  
  
  if (POLARIZATION=='VV') { s1 = s1.select(['VV','angle'])}
  else if (POLARIZATION=='VH') {s1 = s1.select(['VH','angle'])}
  else if (POLARIZATION=='VVVH') {s1 = s1.select(['VV', 'VH', 'angle'])}
  
  print('Number of images in collection: ', s1.size());
  
  /*************** 
   * 2. MASK EDGES  
  ****************/
  
  if (APPLY_BORDER_NOISE_CORRECTION===true) {var s1_1 = s1.map(border_noise_correction.f_mask_edges)}
  if (APPLY_BORDER_NOISE_CORRECTION===false) {s1_1 = s1}


  /*************************  
   * 3. SPECKLE FILTER  
  *************************/
 if (APPLY_SPECKLE_FILTERING) {
    if (SPECKLE_FILTER_FRAMEWORK == 'MONO') {
        s1_1 = ee.ImageCollection(speckle_filter.MonoTemporal_Filter(s1_1, SPECKLE_FILTER_KERNEL_SIZE, SPECKLE_FILTER ))
  }
    else {
        s1_1 = ee.ImageCollection(speckle_filter.MultiTemporal_Filter(s1_1, SPECKLE_FILTER_KERNEL_SIZE, SPECKLE_FILTER,NR_OF_IMAGES ));
  }    
  print('SPECKLE FILTERING COMPLETED') 
 }  
 
   /*************************  
   * 4. TERRAIN CORRECTION  
  **************************/
  
  if (APPLY_TERRAIN_FLATTENING) {
      s1_1 = ee.ImageCollection(terrain_flattening.slope_correction(s1_1, TERRAIN_FLATTENING_MODEL, DEM, TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER)); 
      print('TERRAIN FLATTENING COMPLETED')
  }
  
  return s1_1;
};
