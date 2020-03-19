""" 
RPERefinementAlgorithm: Pipeline to run RPE refinement with larger neighborhood
---------------------------------------------------------------------------------
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
from Algorithms.RefinementRPE import GraphWeightsRefine
from Algorithms.RefinementRPE.GraphCutSegmentation import execute_graphcut
import warnings
from FileHandler.ParameterReader import readRPERefinementDict
warnings.filterwarnings("ignore")

class RPERefinementClass:
    """
        Class to refine the RPE segmentation of OCT images.
        
        Different to the three layer segmentation, this algorithm uses 
        adapted weights and a 9-neighborhood. In addition, the area beneath
        Bruch's membrane is excluded to exclude CC area and area above ILM
        to speed up the process. 
        
        Parameter used from Parameters text:
            - REF_RPE_BF_ENFACE: En-Face smoothing bilateral filter values
            - REF_RPE_BF_BSCAN : B-scan smoothing bilateral filter values
    """
    def __init__(self, volume, segmentation, shiftedValues, scaling = 1):
        """
            Initializing
            
            Parameters
            -----------
            volume: ndarray 
                volume
            segmentation: ndarray 
                segmentation
            Optional
            -----------
            scaling: scalar, float
                up/downscaling factor 
            
        """
        self.scaling = scaling
        max_value = np.max(volume)
        vol = volume.copy()
        #get parameters from parameters text
        self.dictParameters = readRPERefinementDict()
        
        im_center_y = int(segmentation.shape[1]//2)  
        self.shiftedValues = shiftedValues.astype('int32')
        #flattening to Bruch's membrane
        for z in range(vol.shape[0]):
            coordinates_bruchs = np.where(segmentation[z].transpose() == self.dictParameters['BM_VALUE'])[1]
            self.shiftedValues[z,:]= -coordinates_bruchs[:] + im_center_y
            for x in range (vol.shape[2]):
                vol[z,:,x] = np.roll(vol[z,:,x],self.shiftedValues[z][x], axis= 0)  
                   
        for z in range(segmentation.shape[0]):
            for x in range (segmentation.shape[2]):          
                segmentation[z,:,x] = np.roll(segmentation[z,:,x],self.shiftedValues[z][x], axis= 0)
        
        self.oct_volume = [vol_slice/max_value for vol_slice in vol]

        self.segmentation = segmentation
        
    def run(self):
        """
            Execute Multithreaded Pipeline
        """        
        try:
            #Calculate graph weights and run multi-threaded pipeline
            smoothed = GraphWeightsRefine.runWeightCalculation(self.oct_volume, self.segmentation, self.dictParameters)
            self.result = execute_graphcut(self.oct_volume, smoothed, self.scaling, 'RPE', self.dictParameters)
        except Exception as e:
            print('Failed:', e)
        
    def getResult(self):
        """
            Return result of pipeline calculation
            
            Returns
            ---------
            result: ndarray
                inpainted BM, RPE and ILM layer
        """  
        self.result = np.where(self.segmentation == self.dictParameters['BM_VALUE'], self.dictParameters['BM_VALUE'], self.result)
        self.result = np.where(self.segmentation == self.dictParameters['ILM_VALUE'], self.dictParameters['ILM_VALUE'], self.result)
        
        #flatten
        for z in range(self.result.shape[0]):
            for x in range (self.result.shape[2]):          
                self.result[z,:,x] = np.roll(self.result[z,:,x],-self.shiftedValues[z][x], axis= 0)
                
        self.shiftedValues =np.zeros((self.shiftedValues.shape[0], self.shiftedValues.shape[1])).astype('int32')
        return self.result
    
    def runPipeline(self):
        """
            Execute Multithreaded Pipeline
            
            Returns
            ---------
            result: ndarray
                inpainted BM, RPE and ILM layer
        """     
        self.run()
        return self.getResult().astype('uint8')