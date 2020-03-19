""" 
BruchsSegmentation: Module for Approximating Bruch's Membrane
-------------------------------------------------------------
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
from scipy import ndimage
#Fixed Values
MAXITER = 50
BM_OFFSET = 4

def calculate_bruchs(segmentation, dictParameters):
    """ 
        Approximate Bruch's membrane based on RPE
        
        Parameters
        ----------
        segmentation: ndarray
            the underlying segmentation
        dictParameters: dictionary
            Parameters from parameters.txt (here: BM_VALUE and RPE_VALUE)
            
        Returns
        ---------
        result: ndarray
            the resulting segmentation
            
    """
    #get RPE and set up initial BM (RPE+1)
    result_init_coord_rpe_floor = np.where(segmentation == int(dictParameters['RPE_VALUE']))
    result = np.zeros((segmentation.shape))
    result[result_init_coord_rpe_floor[0],result_init_coord_rpe_floor[1],result_init_coord_rpe_floor[2]] = int(dictParameters['RPE_VALUE'])
    result[result_init_coord_rpe_floor[0],result_init_coord_rpe_floor[1] + 1,result_init_coord_rpe_floor[2]] = int(dictParameters['BM_VALUE'])
    
    cur_iter = 0
    eps = np.inf
    #iterate
    while(cur_iter < 10 and eps > 1.0):

        result = bruchs_algorithm(result,dictParameters)
        result = bruchs_algorithm(np.swapaxes(result, axis1= 0 , axis2 = 2),dictParameters)
        result = np.swapaxes(result, axis1= 2 , axis2 = 0)
        eps = np.sqrt(((result - result) ** 2).mean())

        cur_iter += 1
    #write result   
    for z in range (result.shape[0]):
        for x in range (result.shape[2]):
            val = np.where(result[z, : ,x] == int(dictParameters['BM_VALUE']))[0][-1]
            result[z, val ,x] = 0
            result[z, val + BM_OFFSET ,x] = int(dictParameters['BM_VALUE'])
    
    return result
        
def bruchs_algorithm(segmentation, dictParameters):
    """ 
        Helper to approximate Bruch's membrane based on RPE
        
        Parameters
        ----------
        segmentation: ndarray
            the underlying segmentation
        dictParameters: dictionary
            Parameters from parameters.txt (here: BM_VALUE and RPE_VALUE)
          
        Returns
        ---------
        result: ndarray
            the resulting segmentation
            
    """    
    result = np.zeros((segmentation.shape)).astype('uint8')
    bruchs_x = np.arange(segmentation.shape[2]).astype(np.int32).flatten()
    
    #mean filter kernel = 10% of B-scan width
    N = segmentation.shape[2]//10
    box = np.ones(N)/N
    #iterate through volume until convergence: mean filter BM and extract max of BM an RPE+1
    for i,_slice in enumerate(segmentation):
        
        rpe_floor_y = np.zeros(_slice.shape[1]).astype(np.int32)
        rpe_y = np.zeros(_slice.shape[1]).astype(np.int32)
        
        for x in range (_slice.shape[1]):
            rpe_y[x] = np.where(_slice[:,x] == int(dictParameters['RPE_VALUE']))[0][-1]
            rpe_floor_y[x] = np.where(_slice[:,x] == int(dictParameters['BM_VALUE']))[0][-1]

        bruchs_y = rpe_floor_y.copy().flatten()
        rpe_floor_y = rpe_floor_y.flatten()

        eps = np.inf
        cur_iter = 0

        while(cur_iter < MAXITER):

            y = ndimage.convolve(bruchs_y, box, mode='reflect')
            
            bruchs_y_max = np.maximum(y.flatten(), rpe_floor_y.flatten())
            eps = np.sqrt(((bruchs_y_max - bruchs_y) ** 2).mean())
            bruchs_y_max[0:20] = rpe_floor_y[0:20]
            bruchs_y_max[-20:] = rpe_floor_y[-20:]
            bruchs_y = bruchs_y_max.copy()
            
            if(eps < 0.05):
                break
            
            cur_iter += 1
            
        bruchs_y = ndimage.convolve(bruchs_y, box, mode='reflect')
        
        rpe_new_y = np.where(rpe_y >= bruchs_y, bruchs_y - 1, rpe_y) 
        
        result[i,rpe_new_y,bruchs_x] = int(dictParameters['RPE_VALUE'])
        result[i,bruchs_y,bruchs_x] = int(dictParameters['BM_VALUE'])
        
    return result