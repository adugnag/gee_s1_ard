#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar  6 10:55:16 2021

@author: adugnamullissa
"""

import ee
import math
import geetools

ee.Initialize()

#---------------------------------------------------------------------------//
# Linear to db scale
#---------------------------------------------------------------------------//

def lin_to_db(image):
  """ Convert backscatter from linear to dB. """
  bandNames = image.bandNames()
  db = ee.Image(10).multiply(image.log10()).rename(bandNames)
  return db.set('system:time_start', image.get('system:time_start'))


#---------------------------------------------------------------------------//
# Additional Border Noise Removal
#---------------------------------------------------------------------------//


def maskAngLT452(image):
   """/** Function to mask out edges of images using using angle.
   * Source: Hird et al. 2017 Remote Sensing (supplementary material): http://www.mdpi.com/2072-4292/9/12/1315)
   * (mask out angles >= 45.23993) */"""
   ang = image.select(['angle'])
   return image.updateMask(ang.lt(45.23993)).set('system:time_start', image.get('system:time_start'))



def maskAngGT30(image):
   """/** Function to mask out edges of images using angle.
   * (mask out angles <= 30.63993) */"""
   ang = image.select(['angle'])
   return image.updateMask(ang.gt(30.63993)).set('system:time_start', image.get('system:time_start'))


def maskEdge(image):
   """/** Remove edges.
   * Source: Andreas Vollrath */"""
   mask = image.select(0).unitScale(0.000316, 3).multiply(255).toByte().connectedComponents(ee.Kernel.rectangle(1,1), 100)
   return image.updateMask(mask.select(0)).set('system:time_start', image.get('system:time_start')) 



def f_mask_edges(image):
  """/** Mask edges. This function requires that the input image has one VH or VV band, and an 'angle' bands.  */"""
  output = maskAngGT30(image)
  output = maskAngLT452(output)
  output = maskEdge(output)
  return output.set('system:time_start', image.get('system:time_start'))


#---------------------------------------------------------------------------//
# Terrain Flattening
#---------------------------------------------------------------------------//
 
def slope_correction(collection, TERRAIN_FLATTENING_MODEL \
                                              ,DEM \
                                              ,POLARIZATION \
                                              ,TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER):
  
  """ Vollrath, A., Mullissa, A., & Reiche, J. (2020). Angular-Based Radiometric Slope Correction for Sentinel-1 on Google Earth Engine. 
  Remote Sensing, 12(11), [1867]. https://doi.org/10.3390/rs12111867 """
  
  ninetyRad = ee.Image.constant(90).multiply(math.PI/180)

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
      theta_iRad = image.select('angle').multiply(math.PI/180)
      phi_iRad = ee.Image.constant(heading).multiply(math.PI/180)
        
      #2.1.2 Terrain geometry
      alpha_sRad = ee.Terrain.slope(DEM).select('slope').multiply(math.PI/180)#.reproject(proj).clip(geom)
      phi_sRad = ee.Terrain.aspect(DEM).select('aspect').multiply(math.PI/180)#.reproject(proj).clip(geom)
        
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
        
      output = ee.Image(output).addBands(image.select('angle'),'null', True)
        
      return output.set('system:time_start', image.get('system:time_start'))  
  return collection.map(_correct)


#---------------------------------------------------------------------------//
# Mono-temporal speckle filters
#---------------------------------------------------------------------------//

def boxcar(image, KERNEL_SIZE):
  """Applies boxcar filter on every image in the collection."""
  bandNames = image.bandNames()
    #Define a boxcar kernel
  kernel = ee.Kernel.square((KERNEL_SIZE/2), 'pixels')
    #Apply boxcar
  output = image.convolve(kernel)
  output = output.rename(bandNames).copyProperties(image)
  return output.set('system:time_start', image.get('system:time_start'))


def leefilter(image,KERNEL_SIZE):
    """Lee Filter applied to one image. It is implemented as described in 
    J. S. Lee, “Digital image enhancement and noise filtering by use of local statistics,” 
    IEEE Pattern Anal. Machine Intell., vol. PAMI-2, pp. 165–168, Mar. 1980."""
        
    bandNames = image.bandNames()
  
    #S1-GRD images are multilooked 5 times in range
    enl = 5
    #Compute the speckle standard deviation
    eta = 1.0/math.sqrt(enl) 
    eta = ee.Image(eta)

    #MMSE estimator
    #Neighbourhood mean and variance
    oneImg = ee.Image(1)
    z_bar = image.reduceNeighborhood(reducer= ee.Reducer.mean(),kernel= ee.Kernel.square((KERNEL_SIZE/2), 'pixels'))
    varz = image.reduceNeighborhood(reducer= ee.Reducer.variance(),kernel= ee.Kernel.square((KERNEL_SIZE/2), 'pixels'))
    #Estimate weight 
    varx = (varz.subtract(z_bar.pow(2).multiply(eta.pow(2)))).divide(oneImg.add(eta.pow(2)))
    b = varx.divide(varz)
  
    #if b is negative set it to zero
    value = 0 
    #Create mask
    mask = b.lt(0)
    #Create new weight
    new_b = mask.multiply(value).add(b.multiply(mask.Not())) 
    output = oneImg.subtract(new_b).multiply(z_bar.abs()).add(new_b.multiply(image))
    output = output.rename(bandNames).copyProperties(image)
    return output.set('system:time_start', image.get('system:time_start'))


def gammamap(image,KERNEL_SIZE): 
    """ Gamma Maximum a-posterior Filter applied to one image. It is implemented as described in 
    Lopes A., Nezry, E., Touzi, R., and Laur, H., 1990.  Maximum A Posteriori Speckle Filtering and First Order texture Models in SAR Images.  
    International  Geoscience  and  Remote  Sensing  Symposium (IGARSS). """
    enl = 5
    bandNames = image.bandNames()
    #local mean
    z = image.reduceNeighborhood(reducer= ee.Reducer.mean(),kernel= ee.Kernel.square((KERNEL_SIZE/2), 'pixels'))
    #local standard deviation
    sigz = image.reduceNeighborhood(reducer= ee.Reducer.StdDev(),kernel= ee.Kernel.square((KERNEL_SIZE/2), 'pixels'))
    #local observed coefficient of variation
    ci = sigz.divide(z)
    #noise coefficient of variation (or noise sigma)
    cu = 1.0/math.sqrt(enl)
    #threshold for the observed coefficient of variation
    cmax = math.sqrt(2.0) * cu
    cu = ee.Image(cu)
    cmax = ee.Image(cmax)
    enlImg = ee.Image(enl)
    oneImg = ee.Image(1)
    twoImg = ee.Image(2)
    fourImg = ee.Image(4)
    alpha = oneImg.add(cu.pow(2)).divide(ci.pow(2).subtract(cu.pow(2)))

    #Implements the Gamma MAP filter described in equation 11 in Lopez et al. 1990
    q = z.pow(2).multiply(z.multiply(alpha.subtract(enlImg).subtract(oneImg))
            .pow(2)).add(fourImg.multiply(alpha).multiply(enlImg).multiply(image).multiply(z))
    rHat = z.multiply(alpha.subtract(enlImg).subtract(oneImg)).add(q.sqrt()).divide(twoImg.multiply(alpha))
  
    #if ci <= cu then its a homogenous region ->> boxcar filter
    zHat = (z.updateMask(ci.lte(cu))).rename(bandNames)
    #if cmax > ci > cu then its a textured medium ->> apply Gamma MAP filter
    rHat = (rHat.updateMask(ci.gt(cu)).updateMask(ci.lt(cmax))).rename(bandNames)
    #ci>cmax then its strong signal ->> retain
    x = image.updateMask(ci.gte(cmax)).rename(bandNames)  
    #Merge
    output = ee.ImageCollection([zHat,rHat,x]).sum().rename(bandNames).copyProperties(image)
    return output.set('system:time_start', image.get('system:time_start')) 

def RefinedLee(img):
  """ This filter is modified from the implementation by Guido Lemoine 
  Source: Lemoine et al. https://code.earthengine.google.com/5d1ed0a0f0417f098fdfd2fa137c3d0c
  """
  bandNames = img.bandNames()
  # Set up 3x3 kernels 
  weights3 = ee.List.repeat(ee.List.repeat(1,3),3)
  kernel3 = ee.Kernel.fixed(3,3, weights3, 1, 1, False)

  mean3 = img.reduceNeighborhood(ee.Reducer.mean(), kernel3)
  variance3 = img.reduceNeighborhood(ee.Reducer.variance(), kernel3)

  # Use a sample of the 3x3 windows inside a 7x7 windows to determine gradients and directions
  sample_weights = ee.List([[0,0,0,0,0,0,0], [0,1,0,1,0,1,0],[0,0,0,0,0,0,0], [0,1,0,1,0,1,0], [0,0,0,0,0,0,0], [0,1,0,1,0,1,0],[0,0,0,0,0,0,0]])

  sample_kernel = ee.Kernel.fixed(7,7, sample_weights, 3,3, False)

  # Calculate mean and variance for the sampled windows and store as 9 bands
  sample_mean = mean3.neighborhoodToBands(sample_kernel) 
  sample_var = variance3.neighborhoodToBands(sample_kernel)

  # Determine the 4 gradients for the sampled windows
  gradients = sample_mean.select(1).subtract(sample_mean.select(7)).abs()
  gradients = gradients.addBands(sample_mean.select(6).subtract(sample_mean.select(2)).abs())
  gradients = gradients.addBands(sample_mean.select(3).subtract(sample_mean.select(5)).abs())
  gradients = gradients.addBands(sample_mean.select(0).subtract(sample_mean.select(8)).abs())

  # And find the maximum gradient amongst gradient bands
  max_gradient = gradients.reduce(ee.Reducer.max())

  # Create a mask for band pixels that are the maximum gradient
  gradmask = gradients.eq(max_gradient)

  # duplicate gradmask bands: each gradient represents 2 directions
  gradmask = gradmask.addBands(gradmask)

  # Determine the 8 directions
  directions = sample_mean.select(1).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(7))).multiply(1)
  directions = directions.addBands(sample_mean.select(6).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(2))).multiply(2))
  directions = directions.addBands(sample_mean.select(3).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(5))).multiply(3))
  directions = directions.addBands(sample_mean.select(0).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(8))).multiply(4))
  # The next 4 are the not() of the previous 4
  directions = directions.addBands(directions.select(0).Not().multiply(5))
  directions = directions.addBands(directions.select(1).Not().multiply(6))
  directions = directions.addBands(directions.select(2).Not().multiply(7))
  directions = directions.addBands(directions.select(3).Not().multiply(8))

  # Mask all values that are not 1-8
  directions = directions.updateMask(gradmask)

  # "collapse" the stack into a singe band image (due to masking, each pixel has just one value (1-8) in it's directional band, and is otherwise masked)
  directions = directions.reduce(ee.Reducer.sum())  

  #var pal = ['ffffff','ff0000','ffff00', '00ff00', '00ffff', '0000ff', 'ff00ff', '000000']
  #Map.addLayer(directions.reduce(ee.Reducer.sum()), {min:1, max:8, palette: pal}, 'Directions', false)

  sample_stats = sample_var.divide(sample_mean.multiply(sample_mean))

  # Calculate localNoiseVariance
  sigmaV = sample_stats.toArray().arraySort().arraySlice(0,0,5).arrayReduce(ee.Reducer.mean(), [0])

  # Set up the 7*7 kernels for directional statistics
  rect_weights = ee.List.repeat(ee.List.repeat(0,7),3).cat(ee.List.repeat(ee.List.repeat(1,7),4))

  diag_weights = ee.List([[1,0,0,0,0,0,0], [1,1,0,0,0,0,0], [1,1,1,0,0,0,0], [1,1,1,1,0,0,0], [1,1,1,1,1,0,0], [1,1,1,1,1,1,0], [1,1,1,1,1,1,1]])

  rect_kernel = ee.Kernel.fixed(7,7, rect_weights, 3, 3, False)
  diag_kernel = ee.Kernel.fixed(7,7, diag_weights, 3, 3, False)

  # Create stacks for mean and variance using the original kernels. Mask with relevant direction.
  dir_mean = img.reduceNeighborhood(ee.Reducer.mean(), rect_kernel).updateMask(directions.eq(1))
  dir_var = img.reduceNeighborhood(ee.Reducer.variance(), rect_kernel).updateMask(directions.eq(1))

  dir_mean = dir_mean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), diag_kernel).updateMask(directions.eq(2)))
  dir_var = dir_var.addBands(img.reduceNeighborhood(ee.Reducer.variance(), diag_kernel).updateMask(directions.eq(2)))

  # and add the bands for rotated kernels
  for i in range(1, 4):
    dir_mean = dir_mean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), rect_kernel.rotate(i)).updateMask(directions.eq(2*i+1)))
    dir_var = dir_var.addBands(img.reduceNeighborhood(ee.Reducer.variance(), rect_kernel.rotate(i)).updateMask(directions.eq(2*i+1)))
    dir_mean = dir_mean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), diag_kernel.rotate(i)).updateMask(directions.eq(2*i+2)))
    dir_var = dir_var.addBands(img.reduceNeighborhood(ee.Reducer.variance(), diag_kernel.rotate(i)).updateMask(directions.eq(2*i+2)))

  # "collapse" the stack into a single band image (due to masking, each pixel has just one value in it's directional band, and is otherwise masked)
  dir_mean = dir_mean.reduce(ee.Reducer.sum())
  dir_var = dir_var.reduce(ee.Reducer.sum())

  # A finally generate the filtered value
  varX = dir_var.subtract(dir_mean.multiply(dir_mean).multiply(sigmaV)).divide(sigmaV.add(1.0))
  b = varX.divide(dir_var)

  result = dir_mean.add(b.multiply(img.subtract(dir_mean)))
  result = result.arrayFlatten([['sum']]).float().toBands().rename(bandNames).copyProperties(img)
  result = result.set('system:time_start', img.get('system:time_start'))

  return result


def leesigma(image,KERNEL_SIZE):
    """/** Implements the improved lee sigma filter to one image. 
    It is implemented as described in, Lee, J.-S. Wen, J.-H. Ainsworth, T.L. Chen, K.-S. Chen, A.J. Improved sigma filter for speckle filtering of SAR imagery. 
    IEEE Trans. Geosci. Remote Sens. 2009, 47, 202–213. """
    #parameters
    Tk = ee.Image(7) #number of bright pixels in a 3x3 window
    sigma = '0.9'
    enl = '4'
    target_kernel = 3
    bandNames = image.bandNames()
  
    #compute the 98 percentile intensity 
    z98 = ee.Dictionary(image.reduceRegion(
                reducer= ee.Reducer.percentile([98],'null',255,0.001,1e6),
                geometry= image.geometry(),
                scale=10,
                bestEffort='true'
            )).toImage()
  

    #select the strong scatterers to retain
    aboveThresh = image.gte(z98)
    K = aboveThresh.reduceNeighborhood(ee.Reducer.sum()
            ,ee.Kernel.square(target_kernel/2)) 
    retainPixel = K.gte(Tk)
  
  
    #compute the a-priori mean within a 3x3 local window
    #original noise standard deviation
    eta = 1.0/math.sqrt(enl)
    eta = ee.Image(eta)
    #MMSE applied to estimate the apriori mean
    z_bar = image.reduceNeighborhood(reducer= ee.Reducer.mean(),kernel= ee.Kernel.square((target_kernel/2), 'pixels'))
    varz = image.reduceNeighborhood(reducer= ee.Reducer.variance(),kernel= ee.Kernel.square((target_kernel/2), 'pixels'))
    oneImg = ee.Image(1)
    varx = (varz.subtract(z_bar.abs().pow(2).multiply(eta.pow(2)))).divide(oneImg.add(eta.pow(2)))
    b = varx.divide(varz)
    xTilde = oneImg.subtract(b).multiply(z_bar.abs()).add(b.multiply(image))
  
    #step 3: compute the sigma range
    #Lookup table (J.S.Lee et al 2009) for range and eta values for intensity (only 4 look is shown here)
    sigmaLookup = ee.Dictionary({  
            4: ee.Dictionary({
            0.5: ee.Dictionary({'I1': 0.694,'I2': 1.385,'eta': 0.1921}),
            0.6: ee.Dictionary({'I1': 0.630,'I2': 1.495,'eta': 0.2348}),
            0.7: ee.Dictionary({'I1': 0.560,'I2': 1.627,'eta': 0.2825}),
            0.8: ee.Dictionary({'I1': 0.480,'I2': 1.804,'eta': 0.3354}),
            0.9: ee.Dictionary({'I1': 0.378,'I2': 2.094,'eta': 0.3991}),
            0.95: ee.Dictionary({'I1': 0.302,'I2': 2.360,'eta': 0.4391})})})
  
    #extract data from lookup
    looksDict = ee.Dictionary(sigmaLookup.get(ee.String(enl)))
    sigmaImage = ee.Dictionary(looksDict.get(ee.String(sigma))).toImage()
    I1 = sigmaImage.select('I1')
    I2 = sigmaImage.select('I2')
    #new speckle sigma
    nEta = sigmaImage.select('eta')
    #establish the sigma ranges
    I1 = I1.multiply(xTilde)
    I2 = I2.multiply(xTilde)
  
    #step 3: apply MMSE filter for pixels in the sigma range
    #MMSE estimator
    mask = image.gte(I1).Or(image.lte(I2))
    z = image.updateMask(mask)
  
    z_bar = image.reduceNeighborhood(reducer= ee.Reducer.mean(),kernel= ee.Kernel.square((KERNEL_SIZE/2), 'pixels'))
    varz = image.reduceNeighborhood(reducer= ee.Reducer.variance(),kernel= ee.Kernel.square((KERNEL_SIZE/2), 'pixels'))
    varx = (varz.subtract(z_bar.abs().pow(2).multiply(nEta.pow(2)))).divide(oneImg.add(nEta.pow(2)))
    b = varx.divide(varz)
    #if b is negative set it to zero
    value = 0 
    mask1 = b.lt(0)
    #Create new image
    new_b = mask1.multiply(value).add(b.multiply(mask1.Not())) 
    xHat = oneImg.subtract(new_b).multiply(z_bar.abs()).add(new_b.multiply(z))
  
    #remove the applied masks and merge the retained pixels and the filtered pixels
    xHat = image.updateMask(retainPixel).unmask(xHat)
    output = ee.Image(xHat).rename(bandNames).copyProperties(image)
    return output.set('system:time_start', image.get('system:time_start')) 

#------------------------------------//
# Multi-temporal speckle filter helpers
#------------------------------------//

def f_divide(count):
  def inner(image):
    """/*Divides an image by the number of images in the collection */"""
    output = image.divide(count)
    output = output.set('system:time_start', image.get('system:time_start'))
    return output
  
  return inner


def f_multiply(isum):
  def inner(image):
    """/*Multiplies an image by the sum of images in a collection */"""
    output = image.multiply(isum)
    output = output.set('system:time_start', image.get('system:time_start'))
    return output
  return inner

#---------------------------------------------------------------------------//
# 4. Multi-temporal speckle filters
#---------------------------------------------------------------------------//
""" The following Multi-temporal speckle filters are implemented as described in
S. Quegan and J. J. Yu, “Filtering of multichannel SAR images,” 
IEEE Trans Geosci. Remote Sensing, vol. 39, Nov. 2001."""

def MultiTemporal_BoxcarFilter(coll, KERNEL_SIZE):
  """ Multi-temporal boxcar Filter.  """
  def inner(image):
    boxcar_filt = boxcar(image, KERNEL_SIZE)
    output = image.divide(boxcar_filt)
    return output.set('system:time_start', image.get('system:time_start'))
  
  def inner2(image):
      output = boxcar(image, KERNEL_SIZE)
      return output.set('system:time_start', image.get('system:time_start'))

  isum = coll.map(inner).reduce(ee.Reducer.sum())
  count_img = coll.reduce(ee.Reducer.count())
  
  s1_filtered = coll \
                .map(inner2) \
                .map(f_divide(count_img)) \
                .map(f_multiply(isum)) \

  return s1_filtered

def MultiTemporal_LeeFilter(coll, KERNEL_SIZE):
  """ Multi-temporal lee Filter.  """
  def inner(image):
    lee_filt = leefilter(image, KERNEL_SIZE)
    output = image.divide(lee_filt)
    return output.set('system:time_start', image.get('system:time_start'))
  
  def inner2(image):
      output = leefilter(image, KERNEL_SIZE)
      return output.set('system:time_start', image.get('system:time_start'))

  isum = coll.map(inner).reduce(ee.Reducer.sum())
  count_img = coll.reduce(ee.Reducer.count())
  
  s1_filtered = coll \
                .map(inner2) \
                .map(f_divide(count_img)) \
                .map(f_multiply(isum)) \

  return s1_filtered

def MultiTemporal_GammaMAP(coll, KERNEL_SIZE):
  """ Multi-temporal Gamma MAP Filter.  """
  def inner(image):
    gamma_filt = gammamap(image, KERNEL_SIZE)
    output = image.divide(gamma_filt)
    return output.set('system:time_start', image.get('system:time_start'))
  
  def inner2(image):
      output = gammamap(image, KERNEL_SIZE)
      return output.set('system:time_start', image.get('system:time_start'))

  isum = coll.map(inner).reduce(ee.Reducer.sum())
  count_img = coll.reduce(ee.Reducer.count())
  
  s1_filtered = coll \
                .map(inner2) \
                .map(f_divide(count_img)) \
                .map(f_multiply(isum)) \

  return s1_filtered

def MultiTemporal_refinedLee(coll, KERNEL_SIZE):
  """ Multi-temporal refined lee Filter.  """
  def inner(image):
    refinedlee_filt = RefinedLee(image)
    output = image.divide(refinedlee_filt)
    return output.set('system:time_start', image.get('system:time_start'))
  
  def inner2(image):
      output = RefinedLee(image)
      return output.set('system:time_start', image.get('system:time_start'))

  isum = coll.map(inner).reduce(ee.Reducer.sum())
  count_img = coll.reduce(ee.Reducer.count())
  
  s1_filtered = coll \
                .map(inner2) \
                .map(f_divide(count_img)) \
                .map(f_multiply(isum)) \

  return s1_filtered

def MultiTemporal_LeeSigma(coll, KERNEL_SIZE):
  """ Multi-temporal refined lee Filter.  """
  def inner(image):
    refinedlee_filt = leesigma(image,KERNEL_SIZE)
    output = image.divide(refinedlee_filt)
    return output.set('system:time_start', image.get('system:time_start'))
  
  def inner2(image):
      output = leesigma(image, KERNEL_SIZE)
      return output.set('system:time_start', image.get('system:time_start'))

  isum = coll.map(inner).reduce(ee.Reducer.sum())
  count_img = coll.reduce(ee.Reducer.count())
  
  s1_filtered = coll \
                .map(inner2) \
                .map(f_divide(count_img)) \
                .map(f_multiply(isum)) \

  return s1_filtered


#****************************
# DO THE JOB
#****************************

def s1_preproc(params):
    """
    Applies preprocessing to a collection of S1 images. 

    Returns:
        (s1_1) where
            s1_1 (ee.ImageCollection): preprocessed image collection

    """
    
    APPLY_BORDER_NOISE_CORRECTION = params['APPLY_BORDER_NOISE_CORRECTION']
    APPLY_TERRAIN_FLATTENING = params['APPLY_TERRAIN_FLATTENING']
    APPLY_SPECKLE_FILTERING = params['APPLY_SPECKLE_FILTERING']
    POLARIZATION = params['POLARIZATION']
    SPECKLE_FILTER = params['SPECKLE_FILTER']
    SPECKLE_FILTER_KERNEL_SIZE = params['SPECKLE_FILTER_KERNEL_SIZE']
    TERRAIN_FLATTENING_MODEL = params['TERRAIN_FLATTENING_MODEL']
    DEM = params['DEM']
    TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER = params['TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER']
    FORMAT = params['FORMAT']
    START_DATE = params['START_DATE']
    STOP_DATE = params['STOP_DATE']
    ORBIT = params['ORBIT']
    ROI = params['ROI']
    SAVE_ASSET = params['SAVE_ASSET']
    ASSET_ID = params['ASSET_ID']

    #****************************
    # 1. CHECK PARAMETERS
    #****************************
    pol_needed = ['VV', 'VH', 'VVVH']
    if (POLARIZATION not in pol_needed.keys()):
        raise ValueError("ERROR!!! Parameter POLARIZATION not correctly defined")

    model_needed = ['DIRECT', 'VOLUME']
    if (POLARIZATION not in model_needed.keys()):
        raise ValueError("ERROR!!! Parameter TERRAIN_FLATTENING_MODEL not correctly defined")

    format_needed = ['LINEAR', 'DB']
    if (TERRAIN_FLATTENING_MODEL not in format_needed.keys()):
        raise ValueError("ERROR!!! FORMAT not correctly defined")

    format_sfilter = ['MONO BOXCAR', 'MONO LEE', 'MONO GAMMA MAP' \
              ,'MONO REFINED LEE', 'MONO LEE SIGMA', 'MULTI BOXCAR', 'MULTI LEE' \
              ,'MULTI GAMMA MAP', 'MULTI REFINED LEE' \
              ,'MULTI LEE SIGMA']
    if (SPECKLE_FILTER not in format_sfilter.keys()):
        raise ValueError("ERROR!!! SPECKLE_FILTER not correctly defined")

    if (TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER < 0):
        raise ValueError("ERROR!!! TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER not correctly defined")

    if (SPECKLE_FILTER_KERNEL_SIZE <= 0):
        raise ValueError("ERROR!!! SPECKLE_FILTER_KERNEL_SIZE not correctly defined")
    
    
    ###########################################
    # 2. IMPORT COLLECTION
    ###########################################
    
    s1 = ee.ImageCollection('COPERNICUS/S1_GRD_FLOAT') \
                .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH')) \
                .filter(ee.Filter.eq('orbitProperties_pass', ORBIT)) \
                .filter(ee.Filter.eq('instrumentMode', 'IW')) \
                .filter(ee.Filter.eq('resolution_meters', 10)) \
                .filterDate(START_DATE, STOP_DATE) \
                .filterBounds(ROI)
    
    #########################
    # 3. SELECT POLARIZATION
    #########################
    
    if (POLARIZATION == 'VV'):
      s1 = s1.select(['VV','angle'])
    elif (POLARIZATION == 'VH'):
      s1 = s1.select(['VH','angle'])
    elif (POLARIZATION == 'VVVH'):
      s1 = s1.select(['VV','VH','angle'])   
    
    ###########################################
    # 4. ADDITIONAL BORDER NOISE CORRECTION
    ###########################################
    
    if (APPLY_BORDER_NOISE_CORRECTION == True):
      s1_1 = s1.map(f_mask_edges)
    else : 
      s1_1 = s1
      
    ########################
    # 5. TERRAIN CORRECTION
    ####################### 
    
    if (APPLY_TERRAIN_FLATTENING == True):
      s1_1 = slope_correction(s1_1 \
                              , TERRAIN_FLATTENING_MODEL \
                              ,DEM \
                              ,POLARIZATION \
                              ,TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER)
      
      
    ########################
    # 5. SPECKLE FILTERING
    ####################### 
    
    if (APPLY_SPECKLE_FILTERING == True):
        if (SPECKLE_FILTER == 'MONO BOXCAR'):
              def filter(SPECKLE_FILTER_KERNEL_SIZE):
                  def inner(image):
                      output = boxcar(image,SPECKLE_FILTER_KERNEL_SIZE)
                      return output.set('system:time_start', image.get('system:time_start'))
                  return inner
              s1_1 = s1_1.map(filter(SPECKLE_FILTER_KERNEL_SIZE))
        if (SPECKLE_FILTER == 'MONO LEE'):
              def filter(SPECKLE_FILTER_KERNEL_SIZE):
                  def inner(image):
                      output = leefilter(image,SPECKLE_FILTER_KERNEL_SIZE)
                      return output.set('system:time_start', image.get('system:time_start'))
                  return inner
              s1_1 = s1_1.map(filter(SPECKLE_FILTER_KERNEL_SIZE))  
        if (SPECKLE_FILTER == 'MONO GAMMA MAP'):
              def filter(SPECKLE_FILTER_KERNEL_SIZE):
                  def inner(image):
                      output = gammamap(image,SPECKLE_FILTER_KERNEL_SIZE)
                      return output.set('system:time_start', image.get('system:time_start'))
                  return inner
              s1_1 = s1_1.map(filter(SPECKLE_FILTER_KERNEL_SIZE))
        if (SPECKLE_FILTER == 'MONO REFINED LEE'):
              s1_1 = s1_1.map(RefinedLee)
        if (SPECKLE_FILTER == 'MONO LEE SIGMA'):
              def filter(SPECKLE_FILTER_KERNEL_SIZE):
                  def inner(image):
                      output = leesigma(image,SPECKLE_FILTER_KERNEL_SIZE)
                      return output.set('system:time_start', image.get('system:time_start'))
                  return inner
              s1_1 = s1_1.map(filter(SPECKLE_FILTER_KERNEL_SIZE))
        if (SPECKLE_FILTER == 'MULTI BOXCAR'):
              s1_1 = MultiTemporal_BoxcarFilter(s1_1, SPECKLE_FILTER_KERNEL_SIZE)
        if (SPECKLE_FILTER == 'MULTI LEE'):
              s1_1 = MultiTemporal_LeeFilter(s1_1, SPECKLE_FILTER_KERNEL_SIZE)
        if (SPECKLE_FILTER == 'MULTI GAMMA MAP'):
              s1_1 = MultiTemporal_GammaMAP(s1_1, SPECKLE_FILTER_KERNEL_SIZE)
        if (SPECKLE_FILTER == 'MULTI REFINED LEE'):
              s1_1 = MultiTemporal_refinedLee(s1_1, SPECKLE_FILTER_KERNEL_SIZE)
        if (SPECKLE_FILTER == 'MULTI LEE SIGMA'):
              s1_1 = MultiTemporal_LeeSigma(s1_1, SPECKLE_FILTER_KERNEL_SIZE)
              
    ########################
    # 7. FORMAT CONVERSION
    ####################### 
    
    if (FORMAT == 'DB'):
        s1_1 = s1_1.map(lin_to_db)
        
    ########################
    # 8. EXPORT ASSET
    #######################  
    
    if (SAVE_ASSET == True):           

        # Export to Asset
        tasks = geetools.batch.Export.imagecollection.toAsset(
            collection=s1_1,
            assetId=ASSET_ID,
            region=s1_1.first().geometry(),
            scale=10,
            dataType='float',
            maxPixels=1e13
        )
        tasks.start()
        print(tasks.status())
        
    return s1_1
        