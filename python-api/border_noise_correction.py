#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Version: v1.0
Date: 2021-03-12
Authors: Mullissa A., Vollrath A.,  Reiche J., Slagter B., Balling J. , Gou Y., Braun, C.
Description: Function to mask out edges of images using using angle.
    Adopted from: Hird et al. 2017 Remote Sensing (supplementary material): http://www.mdpi.com/2072-4292/9/12/1315)
"""

import ee

#---------------------------------------------------------------------------//
# Additional Border Noise Removal
#---------------------------------------------------------------------------//


def maskAngLT452(image):
   """ (mask out angles >= 45.23993) """
   ang = image.select(['angle'])
   return image.updateMask(ang.lt(45.23993)).set('system:time_start', image.get('system:time_start'))



def maskAngGT30(image):
   """ Function to mask out edges of images using angle.
    (mask out angles <= 30.63993) """
   ang = image.select(['angle'])
   return image.updateMask(ang.gt(30.63993)).set('system:time_start', image.get('system:time_start'))


def maskEdge(image):
   """ Remove edges.
   Source: Andreas Vollrath """
   mask = image.select(0).unitScale(0.000316, 3).multiply(255).toByte().connectedComponents(ee.Kernel.rectangle(1,1), 100)
   return image.updateMask(mask.select(0)).set('system:time_start', image.get('system:time_start')) 



def f_mask_edges(image):
  """ Mask edges. This function requires that the input image has one VH or VV band, and an 'angle' bands. """
  output = maskAngGT30(image)
  output = maskAngLT452(output)
  output = maskEdge(output)
  return output.set('system:time_start', image.get('system:time_start'))
