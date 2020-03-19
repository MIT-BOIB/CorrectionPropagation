""" 
GraphCutSegmentation: Multithreaded Graph-Cut execution
-----------------------------------------------------------
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
import threading
import Algorithms.AutomaticSegmentation.BuildRPEGraph as gm_rpe
import Algorithms.AutomaticSegmentation.BuildILMGraph as gm_vitnfl
from Algorithms.AutomaticSegmentation.InpaintSegmentation import inpaint
#dicts to store values
shortestPath={}
y_offset_dict = {}

class thread_pipeline (threading.Thread):
    """
        Threaded Graph-Cut pipeline
        
        Global Parameters
        ----------------
        shortestPath: dict
            storing the latest shortest path as input for next iteration
        y_offset_dict: dict
            storing the latest y_offset as input for next iteration -> runtime minimization
        
    """
    def __init__(self, threadID ,volume_slice, gradient_slice, prior_shortest_path, scaling, UPORDOWN, mode, dictParameters):
        """
            Initializing Pipeline.
            
            
            Parameters
            ----------
            threadID: scalar
                slice number
            volume_slice: numpy array 2D
                slice that will be processed
            gradient_slice: numpy array 2D
                preprocessed gradient slice
            prior_shortest_path: numpy array 2D
                shortest path of predecessor
            scaling: scalar
                scaling factor of image up or downscaling
            UPORDOWN: string
                direction of algorithm (up/down) - important for storing current shortest path and y offset
            mode: string
                'RPE', or 'ILM' 
            dictParameters: dictionary
                Parameters from parameters.txt
        """
        threading.Thread.__init__(self)
        self.threadID = threadID
        self.volume_slice = volume_slice
        self.gradient_slice = gradient_slice
        self.prior_shortest_path = prior_shortest_path
        self.UPORDOWN = UPORDOWN
        self.scaling = scaling
        self.mode = mode
        self.dictParameters = dictParameters
        
    def run(self):
        """
            Execute Algorithm
        """
        self.runPipeline()
        
    def join(self):
        """
            Returns resulting segmentation slice and thread id
        """
        threading.Thread.join(self)
        return self.segmentation, self.threadID
    
    def runPipeline(self):
        """
            Pipeline implementation
        """
        
        #restore y_offset
        y_offset = y_offset_dict[-1]
        if('down' is self.UPORDOWN):
            y_offset = y_offset_dict[1]
        
        #Calculate graph
        if(self.mode is 'RPE'):
            graph, endnode, y_offset = gm_rpe.getRPEGraph(self.gradient_slice, self.scaling, self.prior_shortest_path, y_offset)
        else:
            graph, endnode, y_offset = gm_vitnfl.getILMGraph(self.gradient_slice, self.scaling, self.prior_shortest_path, y_offset)
            
        #calculate shortest path - if invalid, take predecessor's
        shortest_path = np.asarray(graph.get_shortest_paths(v=0, to=endnode, weights = 'weight'))
        if(shortest_path.size == 0):
            shortest_path = self.prior_shortest_path
        
        #Inpaint segmentation  
        self.segmentation = inpaint(self.volume_slice,shortest_path, y_offset,self.mode, self.dictParameters)

        #store variables for next iteration
        if('down' is self.UPORDOWN):
            shortestPath[1] = shortest_path.flatten()
            y_offset_dict[1] = y_offset
        else:
            shortestPath[-1] = shortest_path.flatten()
            y_offset_dict[-1] = y_offset

    
def execute_graphcut(oct_volume, gradient_volume, scaling, mode, dictParameters):
    """
        Execute Graph-Cut algorithm.
        
        Parameters
        ----------
        oct_volume: numpy array 2D/3D
            list of original volume slices for inpainting the results
        gradient_volume: numpy array 2D/3D
            preprocessed gradient_volume 
        scaling: scalar, optional
            scale factor for image down or upscaling
        mode = string
            mode ='RPE' segments RPE, 'ILM' the inner limiting membrane
        dictParameters: dictionary
            Parameters from parameters.txt
            
        Return
        ------
        result: numpy array 2D/3D
            resulting segmented volume with rpe inpainted as value 127
        
    """
    
    #result array
    result = np.zeros((len(oct_volume),oct_volume[0].shape[0], oct_volume[0].shape[1])).astype("uint8")

    index = 0
    center_slice = int(len(gradient_volume)/2)
    
    #iterate over volume starting from central slice
    while center_slice-index >= 0 or center_slice+index < len(gradient_volume):

        #process initial slice (central)
        if index == 0:
            #compute graph and shortest path
            if(mode is 'RPE'):
                initial_graph, endnode, y_offset = gm_rpe.getRPEGraph(gradient_volume[center_slice], scaling, None, None)
            else:
                initial_graph, endnode, y_offset = gm_vitnfl.getILMGraph(gradient_volume[center_slice], scaling, None, None)
                
            shortestPath[1] = np.asarray(initial_graph.get_shortest_paths(v=0, to=endnode, weights = 'weight')).flatten()  
            shortestPath[-1] = shortestPath[1]
            
            if(shortestPath[1].size == 0):
                raise Exception('no path found')
            #store offset
            y_offset_dict[1] = y_offset
            y_offset_dict[-1] = y_offset

            #inpaint shortest path
            result[center_slice] = inpaint(oct_volume[center_slice],shortestPath[1], y_offset, mode, dictParameters)

            index += 1
            continue

        #multithreaded pipeline execution
        threads = []
        if (center_slice-index) >= 0:
            thread = thread_pipeline((center_slice-index),oct_volume[(center_slice-index)], gradient_volume[(center_slice-index)], shortestPath[1],scaling, 'down', mode, dictParameters)
            thread.start()
            threads.append(thread)
        if (center_slice+index) < len(gradient_volume):
            thread = thread_pipeline((center_slice+index),oct_volume[(center_slice+index)], gradient_volume[(center_slice+index)], shortestPath[-1],scaling, 'up', mode, dictParameters)
            thread.start()
            threads.append(thread)
        
        # get processed slices and wait until all processes are finished
        for t in threads:
            result_part, t_id = t.join()
            result[t_id] = result_part
            
        index = index + 1
        
        # stop when finished
        if(center_slice+index >= len(gradient_volume) and center_slice+index < 0):
            break

    return result