/*
Version: v1.0
Date: 2021-03-12
Authors: Mullissa A., Vollrath A.,  Reiche J., Slagter B., Balling J. , Gou Y., Braun, C.
*/
//****************************
// ALL PREPROCESSING
//****************************

var speckle_filter = require('users/adugnagirma/s1_ard_master:speckle_filter');
var terrain_flattening = require('users/adugnagirma/s1_ard_master:terrain_flattening');
var border_noise_correction = require('users/adugnagirma/s1_ard_master:border_noise_correction');

exports.s1_preproc = function(s1, params) {
  
  /************************  
   * 0. CHECK PARAMETERS  
  ************************/
  if (params.ORBIT === undefined) params.ORBIT = 'BOTH';
  if (params.SPECKLE_FILTER === undefined) params.SPECKLE_FILTER = "GAMMA MAP";
  if (params.SPECKLE_FILTER_KERNEL_SIZE===undefined) params.SPECKLE_FILTER_KERNEL_SIZE = 7;
  if (params.TERRAIN_FLATTENING_MODEL === undefined) params.TERRAIN_FLATTENING_MODEL = 'VOLUME';
  if (params.TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER === undefined) params.TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER = 0;
  if (params.FORMAT === undefined) params.FORMAT = 'DB';
  if (params.POLARIZATION === undefined) params.POLARIZATION = 'VVVH';
  if (params.APPLY_BORDER_NOISE_CORRECTION === undefined) params.APPLY_BORDER_NOISE_CORRECTION = true;
  if (params.APPLY_TERRAIN_FLATTENING===undefined) params.APPLY_TERRAIN_FLATTTENING = true;
  if (params.APPLY_SPECKLE_FILTERING===undefined) params.APPLY_SPECKLE_FILTERING = true; 
  if (params.SPECKLE_FILTER_FRAMEWORK===undefined) params.SPECKLE_FILTER_FRAMEWORK = 'MULTI';
  
  function notContains(list, value) {return list.indexOf(value) == -1;}
  
  var orbit_required = ['ASCENDING', 'DESCENDING', 'BOTH']
  if (notContains(orbit_required, params.ORBIT)) {
       throw new Error("Parameter ORBIT not correctly defined")
  } 
  
  var pol_required = ['VV', 'VH', 'HH', 'HV', 'VVVH', 'HHHV']
  if (notContains(pol_required, params.POLARIZATION)) {
       throw new Error("Parameter POLARIZATION not correctly defined")
  } 
  
  var model_required = ['DIRECT', 'VOLUME']
  if (notContains(model_required, params.TERRAIN_FLATTENING_MODEL)) {
       throw new Error("Parameter TERRAIN_FLATTENING_MODEL not correctly defined")
  } 
  
  var format_required = ['LINEAR', 'DB']
  if (notContains(format_required, params.FORMAT)) {
       throw new Error("Parameter FORMAT not correctly defined")
  } 
  
  var frame_required = ['MONO', 'MULTI']
  if (notContains(frame_required, params.SPECKLE_FILTER_FRAMEWORK)) {
       throw new Error("Parameter SPECKLE_FILTER_FRAMEWORK not correctly defined")
  } 
  
  var sfilter_required = ['BOXCAR', 'LEE', 'GAMMA MAP'
              , 'REFINED LEE', 'LEE SIGMA']
  if (notContains(sfilter_required, params.SPECKLE_FILTER)) {
       throw new Error("Parameter SPECKLE_FILTER not correctly defined")
  } 

  if (params.TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER < 0) {
  throw new Error("The TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER should be GREATER or EQUAL to 0")}

  if (params.SPECKLE_FILTER_KERNEL_SIZE <= 0) {
  throw new Error("The SPECKLE_FILTER_KERNEL_SIZE should be a positive integer")}
  
  /************************  
   * 1. SELECT DATASET
  ************************/ 

  //select orbit
  if (params.ORBIT !== 'BOTH'){var s1 = s1.filter(ee.Filter.eq('orbitProperties_pass', params.ORBIT))}
  
  //select polarization
  if (params.POLARIZATION=='VV') { s1 = s1.select(['VV','angle'])}
  else if (params.POLARIZATION=='VH') {s1 = s1.select(['VH','angle'])}
  else if (params.POLARIZATION=='HH') {s1 = s1.select(['HH','angle'])}
  else if (params.POLARIZATION=='HV') {s1 = s1.select(['HV','angle'])}
  else if (params.POLARIZATION=='VVVH') {s1 = s1.select(['VV', 'VH', 'angle'])}
  else if (params.POLARIZATION=='HHHV') {s1 = s1.select(['HH', 'HV', 'angle'])}
  
  print('Number of images in collection: ', s1.size());
  
  //Clip to roi
  if (params.CLIP_TO_ROI) {s1 = s1.map(function(image) {
              return image.clip(geometry)})}
  
  /************************************  
   * 2. ADDITIONAL BORDER NOISE REMOVAL  
  ****************************** ******/
  
  if (params.APPLY_BORDER_NOISE_CORRECTION) {var s1_1 = s1.map(border_noise_correction.f_mask_edges)}
  else {s1_1 = s1}


  /*************************  
   * 3. SPECKLE FILTER  
  *************************/
 if (params.APPLY_SPECKLE_FILTERING) {
    if (params.SPECKLE_FILTER_FRAMEWORK == 'MONO') {
        s1_1 = ee.ImageCollection(speckle_filter.MonoTemporal_Filter(s1_1, params.SPECKLE_FILTER_KERNEL_SIZE, params.SPECKLE_FILTER ))
  }
    else {
        s1_1 = ee.ImageCollection(speckle_filter.MultiTemporal_Filter(s1_1, params.SPECKLE_FILTER_KERNEL_SIZE, params.SPECKLE_FILTER,params.NR_OF_IMAGES ));
  }    
  print('SPECKLE FILTERING COMPLETED') 
 }  
 
   /***************************************   
   * 4. RADIOMETRIC TERRAIN NORMALIZATION  
  ****************************************/
  
  if (params.APPLY_TERRAIN_FLATTENING) {
      s1_1 = ee.ImageCollection(terrain_flattening.slope_correction(s1_1, params.TERRAIN_FLATTENING_MODEL, params.DEM, params.TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER)); 
      print('TERRAIN FLATTENING COMPLETED')
  }
  
  return s1_1;
};

