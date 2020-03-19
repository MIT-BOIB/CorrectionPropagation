""" 
GraphWeightsRefine: Module to calculate graph weights for RPE refinement
---------------------------------------------------------------------------
PRLEC Framework for OCT Processing and Visualization 
"""
# This framework evolved from a collaboration of:
# - Research Laboratory of Electronics, Massachusetts Institute of Technology, Cambdrige, MA, US
# - Pattern Recognition Lab, Friedrich-Alexander-Universitaet Erlangen-Nuernberg, Germany
# - Department of Biomedical Engineering, Peking University, Beijing, China
# - New England Eye Center, Tufts Medical Center, Boston, MA, US
# v1.0: Updated on Mar 20, 2019
# @author: Daniel Stromer - EMAIL:daniel.stromer@fau.de
# Copyright (C) 2018-2019 - Daniel Stromer
# PRLE is developed as an Open Source project under the GNU General Public License (GPL) v3.0.
import numpy as np
import cv2

def runWeightCalculation(volume, segmentation, dictParameters):
    """
        Calculating gradients for all slices of the volume
        
        The volume is bilateral filtered in the en-face plane. After smoothing each slice with BF, 
        each B-scan is filtered again and area above the ILm and CC area beneath Bruch's Membrane 
        is excluded.
        
        Parameters
        ----------
        volume: ndarray
            oct input volume
        segmentation: ndarray 
            segmentation volume
        
        Return
        ------
        volume_res: list of ndarrays
            list of volume slices that were processed
        
    """
    
    bf_enface = dictParameters['REF_RPE_BF_ENFACE']

    volume_smoothed = np.swapaxes(volume, axis1 = 1, axis2 = 0)
    volume_smoothed = np.array([cv2.bilateralFilter(_slice,bf_enface[0],bf_enface[1],bf_enface[2]) for _slice in volume_smoothed])
    volume_smoothed = np.swapaxes(volume_smoothed, axis1 = 0, axis2 = 1)
    
    volume_res = [calculateGradients(_slice, i, segmentation,dictParameters) for i,_slice in enumerate(volume_smoothed)]
    
    return volume_res

def calculateGradients(_slice, i, segmentation,dictParameters):
    """
        Helper to calculate gradients
        
        Parameters
        ----------
        slice: ndarray 
            oct input volume slice
        i: scalar
            current slice
        segmentation: ndarray
            segmentation volume
        
        Return
        ------
        result: numpy array 2D
            resulting gradient image
    """
    bf_bscan = dictParameters['REF_RPE_BF_BSCAN']
    result = cv2.bilateralFilter(_slice,bf_bscan[0],bf_bscan[1],bf_bscan[2])

    y_coord_BM = np.where(segmentation[i].transpose() == dictParameters['BM_VALUE'])[1]
    y_coord_ILM = np.where(segmentation[i].transpose() == dictParameters['ILM_VALUE'])[1]

    for x in range(segmentation[i].shape[1]):
        result[0:y_coord_ILM[x]+1,x] = -1
        result[y_coord_BM[x]:segmentation[i].shape[1]+1,x] = -1

    return result