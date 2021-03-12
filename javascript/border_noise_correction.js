

//---------------------------------------------------------------------------//
// Additional Border Noise Removal
//---------------------------------------------------------------------------//
/** Function to mask out edges of images using using angle.
 * Source: Hird et al. 2017 Remote Sensing (supplementary material): http://www.mdpi.com/2072-4292/9/12/1315)
 * (mask out angles >= 45.23993) */
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
   
   /** Remove edges.
    * Source: Andreas Vollrath */
   var maskEdge = function(image) {
     var mask = image.select(0).unitScale(0.000316, 3).multiply(255).toByte().connectedComponents(ee.Kernel.rectangle(1,1), 100);
     return image.updateMask(mask.select(0)).set('system:time_start', image.get('system:time_start'));  
   };
   
   /** Mask edges. This function requires that the input image has one VH or VV band, and an 'angle' bands.  */
   exports.f_mask_edges = function(image) {
     var output = maskAngGT30(image);
     output = maskAngLT452(output);
     output = maskEdge(output);
     return output.set('system:time_start', image.get('system:time_start'));
   };
   