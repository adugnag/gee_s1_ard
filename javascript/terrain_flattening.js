
//---------------------------------------------------------------------------//
// Terrain Flattening
//---------------------------------------------------------------------------//
/* Vollrath, A., Mullissa, A., & Reiche, J. (2020). Angular-Based Radiometric Slope Correction for Sentinel-1 on Google Earth Engine. 
Remote Sensing, 12(11), [1867]. https://doi.org/10.3390/rs12111867
*/ 
exports.slope_correction = function(collection, TERRAIN_FLATTENING_MODEL, 
                                              DEM,
                                              TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER) {

  var ninetyRad = ee.Image.constant(90).multiply(Math.PI/180);

  var _volumetric_model_SCF = function(theta_iRad, alpha_rRad) {
      // Volume model
      var nominator = (ninetyRad.subtract(theta_iRad).add(alpha_rRad)).tan();
      var denominator = (ninetyRad.subtract(theta_iRad)).tan();
      return nominator.divide(denominator);
  }

  var _direct_model_SCF = function(theta_iRad, alpha_rRad, alpha_azRad) {
      // Surface model
      var nominator = (ninetyRad.subtract(theta_iRad)).cos();
      var denominator = alpha_azRad.cos()
        .multiply((ninetyRad.subtract(theta_iRad).add(alpha_rRad)).cos());
      return nominator.divide(denominator);
  }
  
  var _erode = function(image, distance)  {
     // buffer function (thanks Noel)
      var d = (image.not().unmask(1)
          .fastDistanceTransform(30).sqrt()
          .multiply(ee.Image.pixelArea().sqrt()));

      return image.updateMask(d.gt(distance));
    }
  
  var _masking = function(alpha_rRad, theta_iRad, buffer){
        // calculate masks
        // layover, where slope > radar viewing angle
        var layover = alpha_rRad.lt(theta_iRad).rename('layover');
        // shadow
        var shadow = alpha_rRad.gt(ee.Image.constant(-1).multiply(ninetyRad.subtract(theta_iRad))).rename('shadow');
        // combine layover and shadow
        var mask = layover.and(shadow);
        // add buffer to final mask
        if (buffer > 0)
            mask = _erode(mask, buffer);
        return mask.rename('no_data_mask');
   }

  var _correct = function(image) {
        var bandNames = image.bandNames();
        // get the image geometry and projection
        var geom = image.geometry()
        var proj = image.select(1).projection()
        
        var elevation = DEM.resample('bilinear').reproject({crs:proj, scale:10}).clip(geom)
        
        // calculate the look direction
        var heading = (ee.Terrain.aspect(image.select('angle'))
                                     .reduceRegion(ee.Reducer.mean(),image.geometry(),1000))
        
        // in case of null values for heading replace with 0
        heading = ee.Dictionary(heading).combine({aspect: 0}, false).get('aspect')
    
        heading = ee.Algorithms.If(
            ee.Number(heading).gt(180),
            ee.Number(heading).subtract(360),
            ee.Number(heading)
        )
        // the numbering follows the article chapters
        // 2.1.1 Radar geometry 
        var theta_iRad = image.select('angle').multiply(Math.PI/180)
        var phi_iRad = ee.Image.constant(heading).multiply(Math.PI/180)
        
        // 2.1.2 Terrain geometry
        //slope 
        var alpha_sRad = ee.Terrain.slope(elevation).select('slope').multiply(Math.PI / 180)

        // aspect (-180 to 180)
        var aspect = ee.Terrain.aspect(elevation).select('aspect').clip(geom)

        // we need to subtract 360 degree from all values above 180 degree
        var aspect_minus = aspect
          .updateMask(aspect.gt(180))
          .subtract(360)

        // we fill the aspect layer with the subtracted values from aspect_minus
        var phi_sRad = aspect
          .updateMask(aspect.lte(180))
          .unmask() 
          .add(aspect_minus.unmask()) //add the minus values
          .multiply(-1)   // make aspect uphill
          .multiply(Math.PI / 180) // make it rad
        
        // we get the height, for export 
        var height = DEM.reproject(proj).clip(geom)
        
        
        // 2.1.3 Model geometry
        //reduce to 3 angle
        var phi_rRad = phi_iRad.subtract(phi_sRad)

        // slope steepness in range (eq. 2)
        var alpha_rRad = (alpha_sRad.tan().multiply(phi_rRad.cos())).atan()

        // slope steepness in azimuth (eq 3)
        var alpha_azRad = (alpha_sRad.tan().multiply(phi_rRad.sin())).atan()

        // local incidence angle (eq. 4)
        var theta_liaRad = (alpha_azRad.cos().multiply((theta_iRad.subtract(alpha_rRad)).cos())).acos()
        var theta_liaDeg = theta_liaRad.multiply(180/Math.PI)

        // 2.2 
        // Gamma_nought
        var gamma0 = image.divide(theta_iRad.cos())

        if (TERRAIN_FLATTENING_MODEL == 'VOLUME') {
            // Volumetric Model
            var scf = _volumetric_model_SCF(theta_iRad, alpha_rRad)
        }
        
        if (TERRAIN_FLATTENING_MODEL == 'DIRECT') {
            var scf = _direct_model_SCF(theta_iRad, alpha_rRad, alpha_azRad)
        }
        // apply model for Gamm0
        var gamma0_flat = gamma0.multiply(scf)

        // get Layover/Shadow mask
        var mask = _masking(alpha_rRad, theta_iRad, TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER);

        var output = gamma0_flat.mask(mask).rename(bandNames).copyProperties(image);
        
        output = ee.Image(output).addBands(image.select('angle'),null,true);
        
        return output.set('system:time_start', image.get('system:time_start')); 
  }   
  return collection.map(_correct)
}
