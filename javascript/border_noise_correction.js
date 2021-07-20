/* File: border_noise_correction.js
Version: v1.1
Date: 2021-03-11
Authors: Adopted from Hird et al. 2017 Remote Sensing (supplementary material): http://www.mdpi.com/2072-4292/9/12/1315)
Description: This script applied additional border noise correction */

var helper = require('users/adugnagirma/gee_s1_ard:utilities');

//---------------------------------------------------------------------------//
// Additional Border Noise Removal
//---------------------------------------------------------------------------//
/** (mask out angles >= 45.23993) */
var maskAngLT452 = function(image) {
 var ang = image.select(['angle']);
 return image.updateMask(ang.lt(45.23993)).set('system:time_start', image.get('system:time_start'));
};

/** Function to mask out edges of images using angle.
 * (mask out angles <= 30.63993) */
var maskAngGT30 = function(image) {
 var ang = image.select(['angle']);
 return image.updateMask(ang.gt(30.63993)).set('system:time_start', image.get('system:time_start'));
};

/** Remove edges.*/
var maskEdge = function(image) {
  var mask = image.select(0).unitScale(-25, 5).multiply(255).toByte()//.connectedComponents(ee.Kernel.rectangle(1,1), 100);
  return image.updateMask(mask.select(0)).set('system:time_start', image.get('system:time_start'));  
};

/** Mask edges. This function requires that the input image has one VH or VV band, and an 'angle' bands.  */
exports.f_mask_edges = function(image) {
  var db_img = helper.lin_to_db(image)
  var output = maskAngGT30(db_img);
  output = maskAngLT452(output);
  //output = maskEdge(output);
  output = helper.db_to_lin(output);
  return output.set('system:time_start', image.get('system:time_start'));
};

