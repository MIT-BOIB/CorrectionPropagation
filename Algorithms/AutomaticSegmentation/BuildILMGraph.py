""" 
BuildILMGraph: Module for building graph of ILM layer
-----------------------------------------------------
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
from igraph import *
import numpy as np
import cv2

w_min = 1e-5

graph_list = {}
endnodes = {}

def getILMGraph(_slice, scaling, shortest_path , y_offset, isMatFile):
    """ 
        Construct the undirected, weighted ILM graph from the weights.
        
        For the initial slice (center), an extensive calculation is conducted.
        Afterwards, the shortest path of the predecessor is dilated to
        speed up the process. Invalid edges are set to np.nan
        
        Parameters
        ----------
        _slice: ndarray
            weight slice
        scaling: scalar
            scaling factor for calculation
        shortest_path: ndarray
            the shortest path from the predecessor
        y_offset: scalar
            the top row extracted from the predecessor
        isMatFile: boolean
            true if mat file, false if tiff file
            
        Returns
        ---------
        g: graph
            the resulting segmentation
        endpoint: scalar, int
            the endpoint of the graph
        y_offset:
            the top row extracted from the predecessor
            
    """
    
    start__slice = 0
    #initial slice
    if shortest_path is None:
        _slice_center_y = _slice.shape[0]//2
        #restrict area of search
        
        #mat file vs tiff file
        start__slice = int((_slice_center_y-300)*scaling)
        end__slice = int((_slice_center_y-60)*scaling)
        if isMatFile == True:
            start__slice = int((_slice_center_y-200)*scaling)
            end__slice = int((_slice_center_y-20)*scaling)
        _slice_buffer = _slice[start__slice:end__slice,:]
        #add border to left and right
        _slice_buffer = cv2.copyMakeBorder(cv2.copyMakeBorder(_slice_buffer, top=0,bottom=0, left=0, right=1, borderType= cv2.BORDER_CONSTANT, value =  2.0), top=0,bottom=0, left=1, right=0, borderType= cv2.BORDER_CONSTANT, value = -2.0)
        _slice_1D = _slice_buffer.flatten()
    #slice has predecessor
    else:
        _slice_buffer = _slice.copy()
        #add border to left and right
        _slice_buffer = cv2.copyMakeBorder(cv2.copyMakeBorder(_slice_buffer, top=0,bottom=0, left=0, right=1, borderType= cv2.BORDER_CONSTANT, value =  2.0), top=0,bottom=0, left=1, right=0, borderType= cv2.BORDER_CONSTANT, value = -2.0)  
        #get predecessors shortest path
        path_x = (shortest_path % _slice_buffer.shape[1]).astype(np.int32)
        path_y = (shortest_path / _slice_buffer.shape[1]).astype(np.int32) + y_offset
        
        start_x = np.where(path_x == 1)[0][-1]
        end_x = np.where(path_x == _slice_buffer.shape[1] - 3)[0][-1]
    
        min_y = np.amin(path_y[int(start_x):int(end_x)])
        max_y = np.amax(path_y[int(start_x):int(end_x)])
    
        buffer_result = np.zeros((_slice_buffer.shape)).astype(np.float32)
        buffer_result[path_y, path_x] = 1.0
        #dilate predecessors shortest path
        #mat file vs tiff file
        kernel = np.maximum(int(10*scaling),2)
        y_offset = np.maximum(int(10*scaling),2)
        if isMatFile == True:
            kernel = np.maximum(int(15*scaling),2)
            y_offset = np.maximum(int(15*scaling),2)

        buffer_result[:,1:-1] = cv2.dilate(buffer_result[:,1:-1], np.ones((kernel,kernel),np.uint8),iterations = 1)
        
        #set invalid pixels to nan
        buffer_result = np.where(buffer_result == 1.0, _slice_buffer, np.nan)
        buffer_result = np.where(_slice_buffer[min_y-y_offset:max_y+y_offset] == 0, 0, buffer_result[min_y-y_offset:max_y+y_offset])
        # set borders to -2/+2
        buffer_result[:,0] = -2
        buffer_result[:,-1] = 2
        _slice_1D = buffer_result.flatten()

    # graph construction
    g = Graph(directed=True)

    g.add_vertices(_slice_1D.size)
    
    sx = _slice_buffer.shape[1]
    edge_list = []
    weight_list = []
    
    for idx in range(_slice_1D.size):
        if np.isnan(_slice_1D[idx]):
            continue
    
        v_i = idx
        
        #pixel is not at the right boundary of the image -> has right neighbors
        if((v_i+1) % sx != 0 ):
            v_j = v_i +1
            if (v_j  >= 0 and v_j  < _slice_1D.shape[0] and np.isnan(_slice_1D[v_j]) == False):
                edge_list.append((v_i, v_j))
                weight_list.append(calculateWeight(_slice_1D[v_i],_slice_1D[v_j]))
                
            v_j = v_i - sx +1
            if (v_j  >= 0 and v_j  < _slice_1D.shape[0] and np.isnan(_slice_1D[v_j]) == False):
                edge_list.append((v_i, v_j))
                weight_list.append(calculateWeight(_slice_1D[v_i],_slice_1D[v_j]))
            v_j = v_i + sx +1
            if (v_j  >= 0 and v_j  < _slice_1D.shape[0] and np.isnan(_slice_1D[v_j]) == False):
                edge_list.append((v_i, v_j))
                weight_list.append(calculateWeight(_slice_1D[v_i],_slice_1D[v_j]))
            
        #pixel is located at the left/right boundary -> zero weights
        if(_slice_1D[v_i] > 1.0 or _slice_1D[v_i] < -1.0):
            _slice_1D[v_i] = 0.0
            v_j = v_i - sx
            if (v_j >= 0 and v_j < _slice_1D.shape[0] and np.isnan(_slice_1D[v_j]) == False):
                edge_list.append((v_i, v_j))
                weight_list.append(w_min)
            v_j = v_i + sx
            if (v_j >= 0 and v_j < _slice_1D.shape[0] and np.isnan(_slice_1D[v_j]) == False):
                edge_list.append((v_i, v_j))
                weight_list.append(w_min)
            
        else:
            #not a boundary pixel
            v_j = v_i - sx
            if (v_j  >= 0 and v_j  < _slice_1D.shape[0] and np.isnan(_slice_1D[v_j]) == False):
                edge_list.append((v_i, v_j))
                weight_list.append(calculateWeight(_slice_1D[v_i],_slice_1D[v_j]))
                
            v_j = v_i + sx
            if (v_j  >= 0 and v_j  < _slice_1D.shape[0] and np.isnan(_slice_1D[v_j]) == False):
                edge_list.append((v_i, v_j))
                weight_list.append(calculateWeight(_slice_1D[v_i],_slice_1D[v_j]))
                
    g.add_edges(edge_list)
    g.es['weight'] = weight_list
    
    if shortest_path is None:
        return g, _slice_1D.shape[0]-1, start__slice
    else:
        return g, _slice_1D.size-1, min_y - y_offset

def calculateWeight(value_v_i, value_v_j):
    """ 
        Exponential weight calculation
        
        When a pixel's intensity is zero, the weight is set to a constant (exp(6))
        
        Parameters
        ----------
        value_v_i: scalar, float
            first pixel's intensity
        value_v_j: scalar, float
            second's pixel's intensity
            
        Returns
        ----------
        weight: scalar, float
            the edge weight
    """
    if(value_v_j > 1.0 or value_v_j < -1.0):
        return np.exp(2.0 - (value_v_i) + w_min)
    elif(value_v_j == 0) or (value_v_i == 0):
        return np.exp(6.0)
    else:
        return np.exp(2.0 - (value_v_i + value_v_j) + w_min)