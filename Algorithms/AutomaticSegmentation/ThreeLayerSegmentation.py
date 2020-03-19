""" 
ThreeLayerSegmentation: Graph-Cut Pipeline for segmenting retinal layers
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
from skimage import io
import os
from Algorithms.AutomaticSegmentation import GraphWeights
from Algorithms.AutomaticSegmentation import BruchsSegmentation as BMSeg
from Algorithms.Flattening.VolumeFlattening import runFlattening, unFlatten
from Algorithms.AutomaticSegmentation.GraphCutSegmentation import execute_graphcut
from FileHandler.ParameterReader import readThreeLayerDict

class LayerSegmentation:
    """
        Class to segment three retinal layers of OCT images (ILM, RPE, BM).
        
        The algorithm uses a graph-cut algorithm on each slice of the volume to 
        segment the RPE boundary (similar to Chiu's method). The result of the 
        prior slice is used as input for the next slice, such that the resulting 
        surface is smooth and accurate. The exponential weights are derived from 
        the smoothed gradient images. 
    """
    def __init__(self, volume, scaling = 1, flattening_factor = 4, statusText=""):
        """
            Initializing
            
            Parameters
            -----------
            volume: ndarray 
                oct volume
                
            Optional
            -----------
            scaling: scalar, float
                up/downscaling factor 
            flattening_factor: scalar, int
                polynomial degree for flattening
            statusText: string variable StringVar()
                Used to track progress in GUI
            
        """
        self.scaling = scaling
        self.statusText=statusText
        self.dictParameters = readThreeLayerDict()
        #set volume to 0...1
        max_value = np.max(volume)
        self.oct_volume = [vol_slice/max_value for vol_slice in volume]
        #Flatten volume to RPE
        try:
            self.statusText.set("Running 3-Layer segmentation...\nFlattening Volume")
        except: 
            pass
        self.oct_volume, self.slice_shifts = runFlattening(self.oct_volume, flattening_factor)
        
    def run(self, mode):
        """
            Running segmentation
            
            Parameters
            -----------
            mode: string 
                'ILM' or 'RPE'
                
        """
        try:
            #Preprocessing volume for GraphCut: Gradient calculation and smoothing
            try:
                self.statusText.set("Running 3-Layer segmentation...\nCalculating Graph Weights of "+mode)
            except:
                pass
            smoothed = GraphWeights.runWeightCalculation(self.oct_volume, mode, self.dictParameters)

            #execute graph cut
            try:
                self.statusText.set("Running 3-Layer segmentation...\nExecuting GraphCut of "+mode)
            except:
                pass
            result = execute_graphcut(self.oct_volume, smoothed, self.scaling, mode, self.dictParameters)

            if mode is 'RPE':
                #Approximate BM from RPE
                try:
                    self.statusText.set("Running 3-Layer segmentation...\nCalculating Bruch's Membrane")
                except:
                    pass
                self.bm_result = BMSeg.calculate_bruchs(result, self.dictParameters)
            elif mode is 'ILM':
                #set ILM result
                self.ilm_result = result.copy()

        except Exception as e:
            print('Failed', e)
        
    def getResult(self):
        """
            Running segmentation
            
            Returns
            -----------
            result: ndarray 
                final three layer segmentation
                
        """
        return np.maximum(self.bm_result, self.ilm_result)
    
    def runPipeline(self):
        """
            Running pipeline for all layers
            
            Returns
            -----------
            result: ndarray 
                final three layer segmentation
                
        """
        #run ILM segmentation
        self.run('ILM')
        #run RPE and BM segmentation
        self.run('RPE')
        #catch result
        result = self.getResult().astype('uint8')
        #unFlatten result and return it
        result = unFlatten(result,self.slice_shifts)
        return result

'''
STANDALONE VERSION
------------------
The main function can be used to run class without GUI

In "path", a base directory is specified (here 'tmp').
In 'tmp', place sub-directories of the data (e.g., patient1, patient2)

The volume name will be identified (here, merged_structural.tiff) and the segmentation run.

Example:

Processing C:/tmp/patient1/merged_structural.tiff
Processing C:/tmp/patient2/merged_structural.tiff
Processing C:/tmp/patient3/merged_structural.tiff
Processing C:/tmp/patient4/merged_structural.tiff
Processing C:/tmp/patient5/merged_structural.tiff

'''

if __name__ == '__main__':
    
    import warnings
    warnings.filterwarnings("ignore")
    
    # path to file
    path = 'C:\\tmp\\' 
    #Batch processing: 
    arr_of_folders = os.listdir(path)
    print('Files to process:', arr_of_folders)
    
    #flattening polynomial
    flattening_polynomial = 4
    
    # for each folder
    for folder in arr_of_folders:
        filepath = path + folder

        arr = os.listdir(filepath)
        for i in range(len(arr)):

            if 'merged_structural.tiff' not in arr[i]:
                continue

            filename_path = filepath +'\\'+ arr[i]
            
            volume = io.imread(filename_path).astype(np.float32)

            layer_seg = LayerSegmentation(volume, scaling = 1, flattening_factor = 4, statusText=None)
            
            result = layer_seg.runPipeline()

            io.imsave(filepath+'\\automated_result.tif',result)