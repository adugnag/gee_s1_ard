/** File: s1_ard_test.js
Version: v1.2
Date: 2021-03-10
Authors: Mullissa A., Vollrath A.,  Reiche J., Slagter B., Balling J. , Gou Y., Braun, C.
Description: This script creates an analysis ready S1 image collection.

    params:
        geometry: The region to include imagery within. (ee.Geometry); The user has to interactively draw a bounding box within the map window or explicitly give coordinates in params.
        ORBITS:  The orbits to include. (string: ASCENDING, DESCENDING OR BOTH)
        START_DATE: The earliest date to include images for (inclusive).
        END_DATE: The latest date to include images for (exclusive).
        POLARIZATION: The Sentinel-1 image polarization to select for processing.
            'VV' - selects the VV polarization.
            'VH' - selects the VH polarization.
            "VVVH' - selects both the VV and VH polarization for processing.
        APPLY_SPECKLE_FILTERING: true or false options to apply speckle filter
        SPECKLE_FILTER: (Optional)    Type of speckle filter to apply.
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
        NR_OF_IMAGES: is the number of multi-temporal images to use in the multi-temporal filter framework.
        APPLY_BORDER_NOISE_CORRECTION: (Optional) true or false options to apply additional Border noise correction:
        TERRAIN_FLATTENING - (Optional) true or false option to apply Terrain correction based on Vollrath et al. 2020 and Hoekman and Reiche 2015. 
        TERRAIN_FLATTENING_MODEL - which model to use for radiometric terrain flattening (DIRECT, or VOLUME)
        DEM - digital elevation model (DEM) to use (as EE asset)
        TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER - buffer parameter for passive layover/shadow mask in meters
        FORMAT - the output format for the processed collection. this can be 'LINEAR' or 'DB'.
        CLIP_TO_ROI - (Optional) Resizes the processed images to the region of interest.
        SAVE_ASSETS - (Optional) Exports the processed collection to an EE asset.
    Returns:
        An ee.ImageCollection with an analysis ready Sentinel 1 imagery with the following bands:
            VV;
            VH;
            angle.
**/

var wrapper = require('users/adugnagirma/s1_ard_master:wrapper');
var helper = require('users/adugnagirma/s1_ard_master:utilities');


//---------------------------------------------------------------------------//
// DEFINE PARAMETERS
//---------------------------------------------------------------------------//

var params = {//1. Data Selection
              START_DATE: "2019-01-01",
              STOP_DATE: "2019-01-02",
              POLARIZATION:'VVVH',
              ORBIT : 'BOTH',
              DEM: ee.Image('USGS/SRTMGL1_003'),
              geometry: ee.Geometry.Polygon([[[104.80, 11.61],[104.80, 11.36],[105.16, 11.36],[105.16, 11.61]]], null, false),
              //2. Additional Border noise correction
              APPLY_BORDER_NOISE_CORRECTION: false,
              //3.Speckle filter
              APPLY_SPECKLE_FILTERING: true,
              SPECKLE_FILTER_FRAMEWORK: 'MULTI',
              SPECKLE_FILTER: 'LEE SIGMA',
              SPECKLE_FILTER_KERNEL_SIZE: 9,
              NR_OF_IMAGES: 10,
              //4. Radiometric terrain normalization
              APPLY_TERRAIN_FLATTENING: true,
              TERRAIN_FLATTENING_MODEL: 'VOLUME',
              TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER: 0,
              //5. Output
              FORMAT : 'DB',
              CLIP_TO_ROI: false,
              SAVE_ASSETS: false
}

//---------------------------------------------------------------------------//
// DO THE JOB
//---------------------------------------------------------------------------//
// Select S1 GRD ImageCollection
var s1 = ee.ImageCollection('COPERNICUS/S1_GRD_FLOAT')
      .filter(ee.Filter.eq('instrumentMode', 'IW'))
      .filter(ee.Filter.eq('resolution_meters', 10))
      .filterDate(params.START_DATE, params.STOP_DATE)
      .filterBounds(params.geometry);

//Preprocess the S1 collection
var s1_preprocces = wrapper.s1_preproc(s1, params);



//---------------------------------------------------------------------------//
// VISUALIZE
//---------------------------------------------------------------------------//

//Visulaization of the first image in the collection in RGB for VV, VH, images
var visparam = {}
if (params.POLARIZATION=='VVVH'){
     if (params.FORMAT=='DB'){
    var s1_preprocces_view = s1_preprocces.map(helper.add_ratio_lin).map(helper.lin_to_db2);
    var s1_view = s1.map(helper.add_ratio_lin).map(helper.lin_to_db2);
    visparam = {bands:['VV','VH','VVVH_ratio'],min: [-20, -25, 1],max: [0, -5, 15]}
    }
    else {
    var s1_preprocces_view = s1_preprocces.map(helper.add_ratio_lin);
    var s1_view = s1.map(helper.add_ratio_lin);
    visparam = {bands:['VV','VH','VVVH_ratio'], min: [0, 0, 0],max: [0.5, 0.5, 20]}
    }
}
else {
    if (params.FORMAT=='DB') {
    s1_preprocces_view = s1_preprocces.map(helper.lin_to_db);
    s1_view = s1.map(helper.lin_to_db);
    visparam = {bands:[params.POLARIZATION],min: -25,max: 0}   
    }
    else {
    s1_preprocces_view = s1_preprocces;
    s1_view = s1;
    visparam = {bands:[params.POLARIZATION],min: 0,max: 0.2}
    }
}

Map.centerObject(params.geometry, 13);

Map.addLayer(s1_view.first(), visparam, 'First image in the input S1 collection', true);
Map.addLayer(s1_preprocces_view.first(), visparam, 'First image in the processed S1 collection', true);

//---------------------------------------------------------------------------//
// EXPORT
//---------------------------------------------------------------------------//

//Convert format for export
if (params.FORMAT=='DB'){
  s1_preprocces = s1_preprocces.map(helper.lin_to_db);
}

//Save processed collection to asset
if(params.SAVE_ASSETS) {
helper.Download.ImageCollection.toAsset(s1_preprocces, '', 
               {scale: 10, 
               region: s1_preprocces.geometry(),
                type: 'float'})
}
