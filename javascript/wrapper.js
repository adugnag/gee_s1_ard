// File: S1_ARD_TEST/wrapper.js
// Version: v1.2
// Date: 2021-04-01
// Authors: Mullissa A., Vollrath A., Braun, C., Slagter B., Balling J., Gou Y., Gorelick N.,  Reiche J.

//****************************
// ALL PREPROCESSING
//****************************

var speckle_filter = require('users/adugnagirma/gee_s1_ard:speckle_filter');
var terrain_flattening = require('users/adugnagirma/gee_s1_ard:terrain_flattening');
var border_noise_correction = require('users/adugnagirma/gee_s1_ard:border_noise_correction');

exports.s1_preproc = function(params) {
  
  /************************  
   * 0. CHECK PARAMETERS  
  ************************/
  if (params.ORBIT === undefined) params.ORBIT = 'BOTH';
  if (params.SPECKLE_FILTER === undefined) params.SPECKLE_FILTER = "GAMMA MAP";
  if (params.SPECKLE_FILTER_KERNEL_SIZE===undefined) params.SPECKLE_FILTER_KERNEL_SIZE = 7;
  if (params.TERRAIN_FLATTENING_MODEL === undefined) params.TERRAIN_FLATTENING_MODEL = 'VOLUME';
  if (params.TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER === undefined) params.TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER = 0;
  if (params.FORMAT === undefined) params.FORMAT = 'DB';
  if (params.DEM === undefined) params.DEM = ee.Image('USGS/SRTMGL1_003');
  if (params.POLARIZATION === undefined) params.POLARIZATION = 'VVVH';
  if (params.APPLY_ADDITIONAL_BORDER_NOISE_CORRECTION === undefined) params.APPLY_ADDITIONAL_BORDER_NOISE_CORRECTION = true;
  if (params.APPLY_TERRAIN_FLATTENING===undefined) params.APPLY_TERRAIN_FLATTTENING = true;
  if (params.APPLY_SPECKLE_FILTERING===undefined) params.APPLY_SPECKLE_FILTERING = true; 
  if (params.SPECKLE_FILTER_FRAMEWORK===undefined) params.SPECKLE_FILTER_FRAMEWORK = 'MULTI';
  
  function notContains(list, value) {return list.indexOf(value) == -1;}
  
  var orbit_required = ['ASCENDING', 'DESCENDING', 'BOTH']
  if (notContains(orbit_required, params.ORBIT)) {
       throw new Error("Parameter ORBIT not correctly defined")
  } 
  
  var pol_required = ['VV', 'VH', 'VVVH']
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
   * 1. Data Selection
  ************************/ 
  
  // Select S1 GRD ImageCollection
var s1 = ee.ImageCollection('COPERNICUS/S1_GRD_FLOAT')
      .filter(ee.Filter.eq('instrumentMode', 'IW'))
      .filter(ee.Filter.eq('resolution_meters', 10))
      .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
      .filterDate(params.START_DATE, params.STOP_DATE)
      .filterBounds(params.GEOMETRY);
  
  //select orbit
  if (params.ORBIT !== 'BOTH'){s1 = s1.filter(ee.Filter.eq('orbitProperties_pass', params.ORBIT))}
  
  //select polarization
  if (params.POLARIZATION=='VV') { s1 = s1.select(['VV','angle'])}
  else if (params.POLARIZATION=='VH') {s1 = s1.select(['VH','angle'])}
  else if (params.POLARIZATION=='VVVH') {s1 = s1.select(['VV', 'VH', 'angle'])}
  
  print('Number of images in collection: ', s1.size());
  
  /************************************  
   * 2. Additional Border Noise Correction  
  ****************************** ******/
  
  if (params.APPLY_ADDITIONAL_BORDER_NOISE_CORRECTION) {
    var s1_1 = s1.map(border_noise_correction.f_mask_edges) 
    print('ADDITIONAL BORDER NOISE CORRECTION COMPLETED') }
  else {s1_1 = s1}


  /*************************  
   * 3. Speckle Filtering  
  *************************/
 if (params.APPLY_SPECKLE_FILTERING) {
    if (params.SPECKLE_FILTER_FRAMEWORK == 'MONO') {
        s1_1 = ee.ImageCollection(speckle_filter.MonoTemporal_Filter(s1_1, params.SPECKLE_FILTER_KERNEL_SIZE, params.SPECKLE_FILTER ))
        print('MONO-TEMPORAL SPECKLE FILTERING COMPLETED') 
  }
    else {
        s1_1 = ee.ImageCollection(speckle_filter.MultiTemporal_Filter(s1_1, params.SPECKLE_FILTER_KERNEL_SIZE, params.SPECKLE_FILTER,params.SPECKLE_FILTER_NR_OF_IMAGES ));
        print('MULTI-TEMPORAL SPECKLE FILTERING COMPLETED') 
  }    
 }  
 
   /***************************************   
   * 4. Radiometric Terrain Normalization 
  ****************************************/
  
  if (params.APPLY_TERRAIN_FLATTENING) {
      s1_1 = ee.ImageCollection(terrain_flattening.slope_correction(s1_1, params.TERRAIN_FLATTENING_MODEL, params.DEM, params.TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER)); 
      print('RADIOMETRIC TERRAIN NORMALIZATION COMPLETED')
  }
  
      //Clip to roi (input)
  if (params.CLIP_TO_ROI) {s1 = s1.map(function(image) {
              return image.clip(params.GEOMETRY)})}
  
        //Clip to roi (processed)
  if (params.CLIP_TO_ROI) {s1_1 = s1_1.map(function(image) {
              return image.clip(params.GEOMETRY)})}
  
  return [s1, s1_1]
};
