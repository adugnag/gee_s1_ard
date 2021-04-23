/* File: s1_ard_test.js
Version: v1.2
Date: 2021-03-10
Authors: Mullissa A., Vollrath A., Gorelick N.,  Reiche J., Slagter B., Balling J. , Gou Y., Braun, C.
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
        DEM : which digital elevation model (DEM) to use (as ee asset)
        TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER : buffer parameter for layover/shadow mask in meters
        FORMAT : the output format for the processed collection. this can be 'LINEAR' or 'DB'.
        SAVE_ASSETS : (Optional) Exports the processed collection to an asset.
    Returns:
        An ee.ImageCollection with an analysis ready Sentinel 1 imagery with the specified polarization images and angle band.
**/

var wrapper = require('users/adugnagirma/s1_ard_master:wrapper');
var helper = require('users/adugnagirma/s1_ard_master:utilities');


//---------------------------------------------------------------------------//
// DEFINE PARAMETERS
//---------------------------------------------------------------------------//

var params = {//1. Data Selection
              START_DATE: "2019-01-01",
              STOP_DATE: "2019-05-30",
              POLARIZATION:'VVVH',
              ORBIT : 'DESCENDING',
              DEM: ee.Image('USGS/SRTMGL1_003'),
              geometry: ee.Geometry.Polygon([[[112.05, -0.25],[112.05, -0.45],[112.25, -0.45],[112.25, -0.25]]], null, false),
              //2. Additional Border noise correction
              APPLY_BORDER_NOISE_CORRECTION: true,
              //3.Speckle filter
              APPLY_SPECKLE_FILTERING: true,
              SPECKLE_FILTER_FRAMEWORK: 'MULTI',
              SPECKLE_FILTER: 'GAMMA MAP',
              SPECKLE_FILTER_KERNEL_SIZE: 9,
              NR_OF_IMAGES: 10,
              //4. Radiometric terrain normalization
              APPLY_TERRAIN_FLATTENING: true,
              TERRAIN_FLATTENING_MODEL: 'VOLUME',
              TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER: 20,
              //5. Output
              FORMAT : 'DB',
              SAVE_ASSETS: false
}

//---------------------------------------------------------------------------//
// DO THE JOB
//---------------------------------------------------------------------------//

//Preprocess the S1 collection
var s1_preprocces = wrapper.s1_preproc(params);

var s1 = s1_preprocces[0]
s1_preprocces = s1_preprocces[1]

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
    visparam = {bands:['VV','VH','VVVH_ratio'], min: [0.01, 0.0032, 1.25],max: [1, 0.31, 31.62]}
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



Map.centerObject(params.geometry, 12);

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
