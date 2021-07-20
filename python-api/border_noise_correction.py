#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Version: v1.1
Date: 2021-03-11
Authors: Adopted from Hird et al. 2017 Remote Sensing (supplementary material): http://www.mdpi.com/2072-4292/9/12/1315)
Description: This script applied additional border noise correction
"""

import ee
import helper

# ---------------------------------------------------------------------------//
# Additional Border Noise Removal
# ---------------------------------------------------------------------------//


def maskAngLT452(image):
    """
    mask out angles >= 45.23993

    Parameters
    ----------
    image : ee.Image
        image to apply the border noise masking

    Returns
    -------
    ee.Image
        Masked image

    """
    ang = image.select(['angle'])
    return image.updateMask(ang.lt(45.23993)).set('system:time_start', image.get('system:time_start'))



def maskAngGT30(image):
    """
    mask out angles <= 30.63993

    Parameters
    ----------
    image : ee.Image
        image to apply the border noise masking

    Returns
    -------
    ee.Image
        Masked image

    """

    ang = image.select(['angle'])
    return image.updateMask(ang.gt(30.63993)).set('system:time_start', image.get('system:time_start'))


def maskEdge(image):
    """
    Remove edges.

    Parameters
    ----------
    image : ee.Image
        image to apply the border noise masking

    Returns
    -------
    ee.Image
        Masked image

    """

    mask = image.select(0).unitScale(-25, 5).multiply(255).toByte()#.connectedComponents(ee.Kernel.rectangle(1,1), 100)
    return image.updateMask(mask.select(0)).set('system:time_start', image.get('system:time_start')) 



def f_mask_edges(image):
    """
    Function to mask out border noise artefacts

    Parameters
    ----------
    image : ee.Image
        image to apply the border noise correction to

    Returns
    -------
    ee.Image
        Corrected image

    """
    
    db_img = helper.lin_to_db(image)
    output = maskAngGT30(db_img)
    output = maskAngLT452(output)
    #output = maskEdge(output)
    output = helper.db_to_lin(output)
    return output.set('system:time_start', image.get('system:time_start'))
