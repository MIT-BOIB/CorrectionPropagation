""" 
Projections: Module to create projections
---------------------------------------------------
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
from skimage.filters import frangi
from scipy.ndimage import median_filter

def getsubRPESlab(volume, segmentation, thickness, projMode):
    """ 
        Calculating sub-RPE slab based on inputs
        
        Parameters
        ----------
        volume: ndarray
            oct volume
        segmentation: ndarray
            segmentation volume
        thickness: scalar, int
            oct volume
        projMode: string
            selected projection mode: min, max, mean, median
            
        Returns
        ----------
        result: ndarray, int32
            2D image of RPEDC distances
        
    """   
    result = np.zeros((segmentation.shape[0],segmentation.shape[2])).astype('uint16')
    for z in range (segmentation.shape[0]):
        for x in range (segmentation.shape[2]):
            #extracting RPE (=64) and starting from 1 below to project with a thickness t
            rpe = np.where(segmentation[z,:,x] == 64)[0][-1] + 1
            if projMode is 'mean':
                result[z,x] = np.mean(volume[z,rpe:rpe+thickness+1,x])
            elif projMode is 'min':
                result[z,x] = np.min(volume[z,rpe:rpe+thickness+1,x])
            elif projMode is 'max':
                result[z,x] = np.max(volume[z,rpe:rpe+thickness+1,x])
            else:
                result[z,x] = np.median(volume[z,rpe:rpe+thickness+1,x])

    return (result - np.min(result))*(65535.0/(np.max(result)-np.min(result)))

def getsubBMSlab(volume, segmentation, thickness, projMode):
    """ 
        Calculating sub-BM slab based on inputs
        
        Parameters
        ----------
        volume: ndarray
            oct volume
        segmentation: ndarray
            segmentation volume
        thickness: scalar, int
            oct volume
        projMode: string
            selected projection mode: min, max, mean, median
            
        Returns
        ----------
        result: ndarray, int32
            2D image of RPEDC distances
        
    """       
    result = np.zeros((segmentation.shape[0],segmentation.shape[2])).astype('uint16')
    for z in range (segmentation.shape[0]):
        for x in range (segmentation.shape[2]):
            #extracting Bruchs membrane (=255) and starting from 1 below to project with a thickness t
            bruchs = np.where(segmentation[z,:,x] == 255)[0][-1] +1
            if projMode is 'mean':
                result[z,x] = np.mean(volume[z,bruchs:bruchs+thickness+1,x])
            elif projMode is 'min':
                result[z,x] = np.min(volume[z,bruchs:bruchs+thickness+1,x])
            elif projMode is 'max':
                result[z,x] = np.max(volume[z,bruchs:bruchs+thickness+1,x])
            else:
                result[z,x] = np.median(volume[z,bruchs:bruchs+thickness+1,x])            

    return (result - np.min(result))*(65535.0/(np.max(result)-np.min(result)))

def getProjectionILMRPE(volume, segmentation, projMode, vesselness = False):
    """ 
        Calculating sub-RPE slab based on inputs
        The algorithm excludes the actual ILM/RPE voxels and starts beneath/above the layers.
        
        Parameters
        ----------
        volume: ndarray
            oct volume
        segmentation: ndarray
            segmentation volume
        projMode: string
            selected projection mode: min, max, mean, median
        vesselness: bool
            vessel calculation or not
            
        Returns
        ----------
        result: ndarray, int32
            2D image of RPEDC distances
        
    """       
    result = np.zeros((segmentation.shape[0],segmentation.shape[2])).astype('uint16')
  
    for z in range (segmentation.shape[0]):
        buf_rpe = 0
        buf_ilm = 0
        for x in range (segmentation.shape[2]):
            #extract ilm and rpe
            ilm    = np.where(segmentation[z,:,x] == 127)[0]
            rpe = np.where(segmentation[z,:,x] == 64)[0]

            if(rpe.size == 0):
                rpe = buf_rpe
            else:
                rpe = rpe[-1]
                buf_rpe = rpe
            
            if(ilm.size == 0):
                ilm = buf_ilm
            else:
                ilm = ilm[-1]
                buf_ilm= ilm
            try:
                #project - the values are chosen such taht the actual rpe and ilm are excluded from the projection
                if projMode is 'mean':
                    result[z,x] = np.mean(volume[z,ilm + 1:rpe,x])
                elif projMode is 'min':
                    result[z,x] = np.min(volume[z,ilm + 1:rpe,x])
                elif projMode is 'max':
                    result[z,x] = np.max(volume[z,ilm + 1:rpe,x])
                elif projMode is 'median':
                    result[z,x] = np.median(volume[z,ilm + 1:rpe,x])   
            except:
                pass
    if(vesselness == True):
        return result

    return (result - np.min(result))*(65535.0/(np.max(result)-np.min(result)))
            
def getVesselnessImage(volume, segmentation):
    """ 
        Calculate a Vessel probabilty map by applying Frangi's algorithm. 
        
        See: https://link.springer.com/chapter/10.1007/BFb0056195
        
        Steps:
        1) Min projection between ILM and RPE
        2) Median filtering the image
        2) Vesselness filtering with scale 1..10 (step 1), beta1=0.5, beta2=0.005
        
        Parameters
        ----------
        volume: ndarray
            oct volume
        segmentation: ndarray
            segmentation volume
            
        Returns
        ----------
        result: ndarray, int32
            2D image of RPEDC distances
        
    """  
    result = median_filter(getProjectionILMRPE(volume, segmentation,'min', True), size=5, mode='reflect')
    return (frangi(result, scale_range=(1, 10), scale_step=1, beta1=0.5, beta2=0.005)*255).astype('uint8')      