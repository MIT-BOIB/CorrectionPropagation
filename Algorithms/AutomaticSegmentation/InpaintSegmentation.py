""" 
InpaintSegmentation: Module to inpaint segmentations into OCT volumes
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

def inpaint(_slice, shortest_path, y_offset, mode, dictParameters):
    """
    Inpaint RPE into volume. Value 64/127 is used to mark the RPE/ILM. 
    
    Parameters
    ----------
    slice: numpy array 2D
        volume slice for getting shape
    shortest_path: numpy array 1D
        shortest path calculated by dijkstra's algorithm
    y_offset: scalar
        current y_offset - denoting the skipped rows for runtime optimization
    mode: string
        'RPE' or 'ILM'
    dictParameters: dictionary
        Parameters from parameters.txt (here: RPE_VALUE and BM_VALUE)
        
    Return
    ------
    segmentation: numpy array 2D
        inpainted segmentation in uint8 format
    
    """
    #add border columns to image
    segmentation = np.zeros((_slice.shape[0], _slice.shape[1]+2)).astype('uint8')
    #convert shortest path from flattened 1-D to 2-D (x and y coordinate)
    path_x = (shortest_path % segmentation.shape[1]).astype(np.int32)
    path_y = (shortest_path / segmentation.shape[1]).astype(np.int32) + y_offset
    #inpaint segmentation line 
    if(mode is 'RPE'):
        segmentation[path_y, path_x] = dictParameters['RPE_VALUE']
    elif(mode is 'ILM'):
        segmentation[path_y, path_x] = dictParameters['ILM_VALUE']
    return segmentation[:,1:-1]