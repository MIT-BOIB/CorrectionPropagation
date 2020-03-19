""" 
GraphCutRPE: Manual Refinement Graph-Cut Pipeline for RPE
-----------------------------------------------------------------
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
import scipy.signal as sp
from FileHandler.ParameterReader import readManRefParameters
#min value for graph weights
w_min = 1e-5

def propagateRPE(volumeOriginal, segmentation, segmentation_original, saved_slices, rect_correction, mode):
    """ 
        Manual refinement algorithm for RPE correction.
        
        The algorithm takes the corrected slices and connects them with lines of certain
        thicknesses in slow-scan direction. The thickness is based on the mode, detected
        by the propagation call. 
        
        If the spacing between subsequent corrected lines is rather low, the thickness of 
        the connections is 1 resulting in the final output. 
        
        If the spacing is rather high, the thickness is t and in this area, a graph cut
        is executed.
        
        Parameter used from Parameters text:
            - MAN_RPE_BF_BSCAN: B-scan smoothing bilateral filter values
            - MAN_RPE_MEDIAN: Median filtering resulting lines
            - MAN_RPE_THICKNESS: Thickness of connecting lines
        
        Parameters
        ----------
        volumeOriginal: ndarray
            original volume
        segmentation: ndarray
            segmentation
        segmentation_original:ndarray
            original segmentation
        saved_slices: dictionary (scalar, ndarray) = (slice number, slice)
            slices that were manually corrected
        rect_correction[0]: scalar
            startpoint x-axis
        rect_correction[2]: scalar
            startpoint y-axis
        rect_correction[1]: scalar
            endpoint x-axis
        rect_correction[3]:scalar
            endpoint y-axis
        mode: string
            'high' or 'low'
            
        Returns
        ---------
        cropped_segmentation: ndarray
            segmentation result in cropped volume
    """
    # get parameters from parameters text
    dictParameters = readManRefParameters()
    #crop volume and swap axes
    part_volume = volumeOriginal[rect_correction[2]-1:rect_correction[3]+1,:,rect_correction[0]-1:rect_correction[1]+2]
    part_volume = part_volume/np.max(part_volume)
    part_volume = np.swapaxes(part_volume, axis1=0, axis2=2)
    cropped_segmentation = np.zeros((rect_correction[3] - rect_correction[2] + 2, segmentation.shape[1], rect_correction[1] - rect_correction[0] + 3)).astype('uint8')
    
    for key,value in saved_slices.items():
        cropped_segmentation[key-rect_correction[2]+1, :,:] = value[:,rect_correction[0]-1:rect_correction[1]+2].astype('uint8')
    #set BM 
    cropped_segmentation = np.where(cropped_segmentation != dictParameters['RPE_VALUE'], 0, cropped_segmentation)
    cropped_segmentation[:,:, 0]  = np.where(segmentation[rect_correction[2]-1:rect_correction[3]+1,: ,rect_correction[0] - 1] == dictParameters['RPE_VALUE'], dictParameters['RPE_VALUE'], 0)
    cropped_segmentation[:,:, -1] = np.where(segmentation[rect_correction[2]-1:rect_correction[3]+1,:,  rect_correction[1] + 1] == dictParameters['RPE_VALUE'], dictParameters['RPE_VALUE'], 0)
    cropped_segmentation[0 ,:,:]  = np.where(segmentation[rect_correction[2]-1,:,rect_correction[0]-1:rect_correction[1]+2] == dictParameters['RPE_VALUE'], dictParameters['RPE_VALUE'], 0)
    cropped_segmentation[-1,:,:]  = np.where(segmentation[rect_correction[3]+1,:,rect_correction[0]-1:rect_correction[1]+2] == dictParameters['RPE_VALUE'], dictParameters['RPE_VALUE'], 0)
        
    cropped_segmentation = np.swapaxes(cropped_segmentation, axis1=0, axis2=2)
    
    #connect points with lines
    cropped_segmentation = connectpoints(cropped_segmentation,mode,dictParameters)

    #Graph-Cut, if too little number of lines for region
    if(mode is 'high'):
        for i in range(1,cropped_segmentation.shape[0]-1):
                graph, endnode = getGraph(part_volume[i], cropped_segmentation[i],dictParameters)
                shortest_path = np.asarray(graph.get_shortest_paths(v=0, to=endnode, weights = 'weight'))
                cropped_segmentation[i] = inpaint(cropped_segmentation[i],shortest_path,dictParameters)

    cropped_segmentation = np.swapaxes(cropped_segmentation, axis1=2, axis2=0) 
    #median filter
    if(mode is 'high'):
        try:
            for i in range(1,cropped_segmentation.shape[0]-1):
                result = np.zeros(cropped_segmentation[i].shape)
                values = np.where(cropped_segmentation[i].transpose() == dictParameters['RPE_VALUE'])
                y_values = sp.medfilt(values[1],dictParameters['MAN_RPE_MEDIAN']).astype('uint16')
                result[y_values,values[0]] = dictParameters['RPE_VALUE']
                cropped_segmentation[i] = result
        except:
            return cropped_segmentation
    #Inpaint RPE and BM
    vals = np.where(np.logical_and(segmentation_original[rect_correction[2]-1:rect_correction[3]+1,:,rect_correction[0]-1:rect_correction[1]+2] == dictParameters['BM_VALUE'], cropped_segmentation == dictParameters['RPE_VALUE']))
    cropped_segmentation = np.where(segmentation_original[rect_correction[2]-1:rect_correction[3]+1,:,rect_correction[0]-1:rect_correction[1]+2] == dictParameters['BM_VALUE'] , dictParameters['BM_VALUE'], cropped_segmentation)
    for i in range(len(vals[1])):
        cropped_segmentation[vals[0][i],vals[1][i]-1,vals[2][i]] = dictParameters['RPE_VALUE']
        cropped_segmentation[vals[0][i],vals[1][i],vals[2][i]] = dictParameters['BM_VALUE']
    #Inpaint ILM
    cropped_segmentation = np.where(segmentation_original[rect_correction[2]-1:rect_correction[3]+1,:,rect_correction[0]-1:rect_correction[1]+2] == dictParameters['ILM_VALUE'] ,dictParameters['ILM_VALUE'], cropped_segmentation).astype('uint8')
    
    return cropped_segmentation

def connectpoints(cropped_segmentation,mode,dictParameters):
    """ 
        Connecting points by line with certain thickness.
        
        The thickness depends on the mode (high or low).
        
        Parameters
        ----------
        cropped_segmentation: ndarray
            cropped segmentation volume
        mode: string
            'high' or 'low'
        dictParameters: dictionary
            Parameters from parameter.txt
            
        Returns
        ---------
        cropped_segmentation: ndarray
            segmentation result in cropped volume
    """    
    for k in range (1,cropped_segmentation.shape[0]-1):
        points = np.asarray(np.where(cropped_segmentation[k].transpose() == dictParameters['RPE_VALUE'])).astype('int32')
        buffer_slice = cropped_segmentation[k].copy()
        
        for i in range (points[0].shape[0]-1):
            if(np.abs(points[0][i] - points[0][i+1]) > 1):
                if(mode is 'low'):
                    cv2.line(buffer_slice,(points[0][i],points[1][i]),(points[0][i+1],points[1][i+1]),dictParameters['RPE_VALUE'],1)
                else:
                    cv2.line(buffer_slice,(points[0][i],points[1][i]),(points[0][i+1],points[1][i+1]),dictParameters['RPE_VALUE'],dictParameters['MAN_RPE_THICKNESS'])
                buffer_slice[:,points[0][i]] = 0
                buffer_slice[points[1][i],points[0][i]] = dictParameters['RPE_VALUE']
        
        buffer_slice[:,points[0][-1]] = 0
        buffer_slice[points[1][-1],points[0][-1]] = dictParameters['RPE_VALUE']
        cropped_segmentation[k] = buffer_slice
    return cropped_segmentation

def getGraph(_slice, shortest_path,dictParameters):
    """ 
        Construct Graph for Manual BM Refinement.
        
        7-neighborhood.
        
        Parameters
        ----------
        _slice: ndarray
            calcualted weights
        shortest_path: ndarray
            shortest path from predecessor
        dictParameters: dictionary
            Parameters from parameter.txt
            
        Returns
        ---------
        g: graph
            resulting graph
        endpoint: scalar
            endpoint of graph
    """ 
    slice_copy = _slice.copy()
    buffer_result = np.zeros((_slice.shape[0],_slice.shape[1]+2))
    buffer_result[:,0] = -2
    buffer_result[:,-1] = 2
    buffer_result[:,1:-1] = shortest_path
    
    bf_bscan = dictParameters['MAN_RPE_BF_BSCAN']
    
    # Gradient filter with kernel width 3 to enforce smoothness
    bf_filtered = cv2.bilateralFilter(slice_copy,bf_bscan[0],bf_bscan[1],bf_bscan[2])
    # allow only areas where a line was inpainted
    buffer_result[:,1:-1] = np.where(buffer_result[:,1:-1] == dictParameters['RPE_VALUE'], bf_filtered, np.nan)
    buffer_result[:,0] = -2
    buffer_result[:,-1] = 2

    #2-D to 1-D conversion
    slice_1D = buffer_result.flatten()
    
    # set up a graph where every pixel is a node (edges are transitions in between)
    g = Graph(directed=True)
    g.add_vertices(slice_1D.size)
    
    sx = buffer_result.shape[1]
    edge_list = []
    weight_list = []

    for idx in range(slice_1D.size):
        if np.isnan(slice_1D[idx]):
            continue
    
        v_i = idx
        #7-pixel neighborhood, to increase/decrease, change here
        if((v_i+1) % sx != 0 ):
            v_j = v_i +1
            if (v_j  >= 0 and v_j  < slice_1D.shape[0] and np.isnan(slice_1D[v_j]) == False):
                edge_list.append((v_i, v_j))
                weight_list.append(calculateWeightExp(slice_1D[v_i],slice_1D[v_j]))
                
            v_j = v_i - sx +1
            if (v_j  >= 0 and v_j  < slice_1D.shape[0] and np.isnan(slice_1D[v_j]) == False):
                edge_list.append((v_i, v_j))
                weight_list.append(calculateWeightExp(slice_1D[v_i],slice_1D[v_j]))

            v_j = v_i + sx +1
            if (v_j  >= 0 and v_j  < slice_1D.shape[0] and np.isnan(slice_1D[v_j]) == False):
                edge_list.append((v_i, v_j))
                weight_list.append(calculateWeightExp(slice_1D[v_i],slice_1D[v_j]))
                
            v_j = v_i - 2*sx +1
            if (v_j  >= 0 and v_j  < slice_1D.shape[0] and np.isnan(slice_1D[v_j]) == False):
                edge_list.append((v_i, v_j))
                weight_list.append(calculateWeightExp(slice_1D[v_i],slice_1D[v_j]))

            v_j = v_i + 2*sx +1
            if (v_j  >= 0 and v_j  < slice_1D.shape[0] and np.isnan(slice_1D[v_j]) == False):
                edge_list.append((v_i, v_j))
                weight_list.append(calculateWeightExp(slice_1D[v_i],slice_1D[v_j]))
                
            v_j = v_i - 3*sx +1
            if (v_j  >= 0 and v_j  < slice_1D.shape[0] and np.isnan(slice_1D[v_j]) == False):
                edge_list.append((v_i, v_j))
                weight_list.append(calculateWeightExp(slice_1D[v_i],slice_1D[v_j]))

            v_j = v_i + 3*sx +1
            if (v_j  >= 0 and v_j  < slice_1D.shape[0] and np.isnan(slice_1D[v_j]) == False):
                edge_list.append((v_i, v_j))
                weight_list.append(calculateWeightExp(slice_1D[v_i],slice_1D[v_j]))
                
            v_j = v_i - 4*sx +1
            if (v_j  >= 0 and v_j  < slice_1D.shape[0] and np.isnan(slice_1D[v_j]) == False):
                edge_list.append((v_i, v_j))
                weight_list.append(calculateWeightExp(slice_1D[v_i],slice_1D[v_j]))

            v_j = v_i + 4*sx +1
            if (v_j  >= 0 and v_j  < slice_1D.shape[0] and np.isnan(slice_1D[v_j]) == False):
                edge_list.append((v_i, v_j))
                weight_list.append(calculateWeightExp(slice_1D[v_i],slice_1D[v_j]))

        if(slice_1D[v_i] > 1.0 or slice_1D[v_i] < -1.0):
            slice_1D[v_i] = 0.0
            v_j = v_i - sx
            if (v_j >= 0 and v_j < slice_1D.shape[0] and np.isnan(slice_1D[v_j]) == False):
                edge_list.append((v_i, v_j))
                weight_list.append(w_min)
            v_j = v_i + sx
            if (v_j >= 0 and v_j < slice_1D.shape[0] and np.isnan(slice_1D[v_j]) == False):
                edge_list.append((v_i, v_j))
                weight_list.append(w_min)
            
        else:
            v_j = v_i - sx
            if (v_j  >= 0 and v_j  < slice_1D.shape[0] and np.isnan(slice_1D[v_j]) == False):
                edge_list.append((v_i, v_j))
                weight_list.append(calculateWeightExp(slice_1D[v_i],slice_1D[v_j]))
                
            v_j = v_i + sx
            if (v_j  >= 0 and v_j  < slice_1D.shape[0] and np.isnan(slice_1D[v_j]) == False):
                edge_list.append((v_i, v_j))
                weight_list.append(calculateWeightExp(slice_1D[v_i],slice_1D[v_j]))
                
    g.add_edges(edge_list)
    g.es['weight'] = weight_list

    return g, slice_1D.size-1
        
def calculateWeightExp(value_v_i, value_v_j):    
    """ 
        Calculate the linear weight between two vertices vi and vj.
        
        Parameters
        ----------
        value_v_i: scalar, float
            intensity at vertices vi
        value_v_j: scalar, float
            intensity at vertices vj
            
        Returns
        ---------
        weight: scalar, float
            linear weight between edge vi and vj
    """ 
    if(value_v_j > 1.0 or value_v_j < -1.0):
        return np.exp(2.0 - (value_v_i) + w_min)
    else:
        return np.exp(2.0 - (value_v_i + value_v_j) + w_min)
    
def inpaint(_slice, shortest_path,dictParameters):
    """ 
        Inpaint path into segmentation.
        
        Parameters
        ----------
        _slice: ndarray
            calcualted weights
        shortest_path: ndarray
            shortest path from predecessor
        dictParameters: dictionary
            Parameters from parameter.txt
            
        Returns
        ---------
        segmentation: ndarray
            resulting inpainted segmentation
    """ 
    segmentation = np.zeros((_slice.shape[0], _slice.shape[1]+2)).astype('uint8')

    path_x = (shortest_path % segmentation.shape[1]).astype(np.int32)
    path_y = (shortest_path / segmentation.shape[1]).astype(np.int32) 

    segmentation[path_y, path_x] = dictParameters['RPE_VALUE']

    return segmentation[:,1:-1]