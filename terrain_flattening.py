#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 11 19:22:58 2021

@author: adugnamullissa
"""

import ee
import math

#---------------------------------------------------------------------------//
# Terrain Flattening
#---------------------------------------------------------------------------//
 
def slope_correction(collection, TERRAIN_FLATTENING_MODEL \
                                              ,DEM \
                                              ,POLARIZATION \
                                              ,TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER):
  
  """ Vollrath, A., Mullissa, A., & Reiche, J. (2020). Angular-Based Radiometric Slope Correction for Sentinel-1 on Google Earth Engine. 
  Remote Sensing, 12(11), [1867]. https://doi.org/10.3390/rs12111867 """
  
  ninetyRad = ee.Image.constant(90).multiply(math.pi/180)

  def _volumetric_model_SCF(theta_iRad, alpha_rRad):

      # Volume model
      nominator = (ninetyRad.subtract(theta_iRad).add(alpha_rRad)).tan()
      denominator = (ninetyRad.subtract(theta_iRad)).tan()
      return nominator.divide(denominator)


  def _direct_model_SCF(theta_iRad, alpha_rRad, alpha_azRad):
      # Surface model
      nominator = (ninetyRad.subtract(theta_iRad)).cos()
      denominator = alpha_azRad.cos() \
        .multiply((ninetyRad.subtract(theta_iRad).add(alpha_rRad)).cos())
      return nominator.divide(denominator)

  
  def _erode(image, distance):
      #buffer function (thanks Noel)
      
      d = (image.Not().unmask(1)
          .fastDistanceTransform(30).sqrt()
          .multiply(ee.Image.pixelArea().sqrt()))
    
      return image.updateMask(d.gt(distance))
    
  
  def _masking(alpha_rRad, theta_iRad, buffer):
      # calculate masks
      # layover, where slope > radar viewing angle
      layover = alpha_rRad.lt(theta_iRad).rename('layover')
      # shadow
      shadow = alpha_rRad.gt(ee.Image.constant(-1).multiply(ninetyRad.subtract(theta_iRad))).rename('shadow')
      # combine layover and shadow
      mask = layover.And(shadow)
      # add buffer to final mask
      if (buffer > 0):
          mask = _erode(mask, buffer)
      return mask.rename('no_data_mask')

  def _correct(image):
        
      bandNames = image.bandNames()
      #get the image geometry and projection
      #geom = image.geometry()
      #proj = image.select(1).projection()
        
      #calculate the look direction
      heading = (ee.Terrain.aspect(image.select('angle'))
                                     .reduceRegion(ee.Reducer.mean(),image.geometry(),1000)
                                     .get('aspect'))
        
      #the numbering follows the article chapters
      #2.1.1 Radar geometry 
      theta_iRad = image.select('angle').multiply(math.pi/180)
      phi_iRad = ee.Image.constant(heading).multiply(math.pi/180)
        
      #2.1.2 Terrain geometry
      alpha_sRad = ee.Terrain.slope(DEM).select('slope').multiply(math.pi/180)#.reproject(proj).clip(geom)
      phi_sRad = ee.Terrain.aspect(DEM).select('aspect').multiply(math.pi/180)#.reproject(proj).clip(geom)
        
      #2.1.3 Model geometry
      #reduce to 3 angle
      phi_rRad = phi_iRad.subtract(phi_sRad)

      #slope steepness in range (eq. 2)
      alpha_rRad = (alpha_sRad.tan().multiply(phi_rRad.cos())).atan()

      #slope steepness in azimuth (eq 3)
      alpha_azRad = (alpha_sRad.tan().multiply(phi_rRad.sin())).atan()

      #2.2 
      #Gamma_nought
      gamma0 = image.divide(theta_iRad.cos())

      if (TERRAIN_FLATTENING_MODEL == 'VOLUME'):
            #Volumetric Model
            scf = _volumetric_model_SCF(theta_iRad, alpha_rRad)
        
      if (TERRAIN_FLATTENING_MODEL == 'DIRECT'):
            scf = _direct_model_SCF(theta_iRad, alpha_rRad, alpha_azRad)

        
       #apply model for Gamm0
      gamma0_flat = gamma0.divide(scf)

      #get Layover/Shadow mask
      mask = _masking(alpha_rRad, theta_iRad, TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER)

      output = gamma0_flat.mask(mask).rename(bandNames).copyProperties(image)
        
      output = ee.Image(output).addBands(image.select('angle'))
        
      return output.set('system:time_start', image.get('system:time_start'))  
  return collection.map(_correct)


