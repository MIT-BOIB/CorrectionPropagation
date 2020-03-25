""" 
GraphWeights: ILM/RPE Graph-Weight calculation for Three layer Segmentation
-----------------------------------------------------------------------------
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
from skimage import exposure

def runWeightCalculation(volume, mode, dictParameters, isMatFile):
    """
        Calculating weights for modes
        
        The volume is bilateral filter in the en-face plane. 
        Next, each B-scan is bilateral filtered again.
        
        RPE: The result of this is returned as weights
        ILM: The gradient is calculated (Dark-to-bright) and the resulting
        slice exposured to stretch the grayvalues to 0...1.
        
        Parameters from parameters:
            - AUTO_ILM_BF_ENFACE : En face plane smoothing values for ILM
            - AUTO_ILM_BF_BSCAN  : B-scan smoothing values for ILM
            - AUTO_RPE_BF_ENFACE : En face plane smoothing values for RPE
            - AUTO_RPE_BF_BSCAN  : B-scan smoothing values for RPE
            
        Parameters
        ----------
        volume: numpy array 3D
            oct input volume
        mode: string
            'RPE' or 'ILM' mode setter to detect filter values
        dictParameters: dictionary
            Parameters from parameters.txt
        
        Return
        ------
        volume: list of numpy arrays
            list of volume slices that were preprocessed
        
    """
    if mode is 'ILM':
        BF_enface = dictParameters['AUTO_ILM_BF_ENFACE']
        BF_bscan = dictParameters['AUTO_ILM_BF_BSCAN']
    else:
        BF_enface = dictParameters['AUTO_RPE_BF_ENFACE']
        BF_bscan = dictParameters['AUTO_RPE_BF_BSCAN']
    
    volume_smoothed = np.swapaxes(volume, axis1 = 1, axis2 = 0)

    volume_smoothed = np.array([cv2.bilateralFilter(_slice,BF_enface[0],BF_enface[1],BF_enface[2]) for _slice in volume_smoothed])

    volume_smoothed = np.swapaxes(volume_smoothed, axis1 = 0, axis2 = 1)
    
    volume_res = [calculateGradients(_slice, mode, BF_bscan,isMatFile) for _slice in volume_smoothed]
    
    return volume_res

def calculateGradients(_slice, mode, BF_bscan,isMatFile):
    """
        Helper to calculate gradients
        
        Parameters
        ----------
        _slice: numpy array 2D
            oct input volume slice
        BF_bscan: ndarray, int, size 3
            Bilteral filter values for smoothing
        
        Return
        ------
        result: numpy array 2D
            resulting gradient image
    """
    result = cv2.bilateralFilter(_slice,BF_bscan[0],BF_bscan[1],BF_bscan[2])

    if(mode is 'ILM'):
        grad = cv2.filter2D(result,-1, np.array([[-1],[1]]))
        grad = np.where(grad > 0 , grad, 0)
    
        img_center = _slice.shape[0]//2
        #mat file vs MIT file
        if isMatFile:
            p2, p98 = np.percentile(grad[img_center-200:img_center-20,:], (2, 98))
            grad[img_center-200:img_center-20,:] = exposure.rescale_intensity(grad[img_center-200:img_center-20,:], in_range=(p2, p98))  
        else:
            p2, p98 = np.percentile(grad[img_center-130:img_center-20,:], (2, 98))
            grad[img_center-130:img_center-20,:] = exposure.rescale_intensity(grad[img_center-130:img_center-20,:], in_range=(p2, p98))
        return grad
    else:
        return result