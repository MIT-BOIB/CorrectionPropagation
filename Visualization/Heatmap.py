""" 
Heatmap: Module to calculate heatmap
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
from tkinter import *
import tkinter.ttk as ttk
                
def calculateHeatmap(segmentation):
    """ 
        Calculating RPEDC distance map from a given segmentation
        
        Parameters
        ----------
        segmentation: ndarray
            given segmentation as tif file
            
        Returns
        ----------
        heatmap_rpedc: ndarray, int32
            2D image of RPEDC distances
        
    """    
    heatmap_rpedc = np.zeros((segmentation.shape[0],segmentation.shape[2])).astype('int32')
    if(np.count_nonzero(segmentation) != 0):
        for z in range (segmentation.shape[0]):
            for x in range (segmentation.shape[2]):
                
                rpe    = np.where(segmentation[z,:,x] == 64)[0]
                bruchs = np.where(segmentation[z,:,x] == 255)[0]
                
                try:
                    val = bruchs - rpe - 5
                    if val >= 0:
                        heatmap_rpedc[z,x] = val
                    else:
                        heatmap_rpedc[z,x] = 0
                except:
                    heatmap_rpedc[z,x] = 0
                    
    heatmap_rpedc -= np.min(heatmap_rpedc)
    return heatmap_rpedc

def calculateMetrics(heatmap, threshold, voxelsize_transversal,voxelsize_axial):
    """ 
        Calculating drusen volumes, color-code decreasing by volume, and create list of volumes.
        There is a bug when it comes to more than 255 drusen (uint8 limit).
        
        Parameters
        ----------
        heatmap: ndarray
            given segmentation as tif file
        threshold: scalar, int
            threshold from slider
        voxelsize_transversal: scalar, float
            transversal resolution
        voxelsize_axial: scalar, float
            axial resolution
            
        Returns
        ----------
        drusenmap: ndarray, int32
            heatmap with colorcoded drusen
        
    """   
    try:
        __import__("cv2")
    except ImportError:
        return print("This module requires OpenCV.")
    import cv2
    
    #threshold heatmap and find connected components
    heatmap = np.where(heatmap < threshold, 0 , heatmap)  
    _, labels = cv2.connectedComponents(heatmap.astype('uint8'),8)
    
    max_label = np.max(labels)+1
    
    #create drusen list
    _sum = 0
    drusenList = []
    for lbl in range(1,max_label):        
        
        labelmap = np.where(labels == lbl,1,0).astype('int32') * heatmap.astype('int32')
        sum_single = np.sum(labelmap)
        
        _sum += sum_single
        
        drusenList.append((lbl,int(np.round(sum_single*voxelsize_transversal**2*voxelsize_axial,0))))
    # sort list
    drusenList = sorted(drusenList,key=lambda x: x[1], reverse=False)
    drusenmap = np.zeros(heatmap.shape).astype('uint8')
    for lbl in range(1, len(drusenList)+1):
        drusenmap = np.where(labels == drusenList[lbl-1][0], lbl, drusenmap)
        
    #create tree view
    metric_root = Tk()
    metricsTree = ttk.Treeview(metric_root, height=40, columns=('No.', 'Volume'), selectmode="extended")
    metricsTree.heading('#0', text='No.', anchor=CENTER)
    metricsTree.heading('#1', text='Volume [\u03BCm\u00b3]', anchor=CENTER)
    metricsTree.column('#0', stretch=YES, minwidth=100, width=100)
    metricsTree.column('#1', stretch=YES, minwidth=100, width=100)
    #insert drusen entries
    lastIdx = len(drusenList)
    metricsTree.insert("", END, text = 'Total',values= int(np.round(_sum*voxelsize_transversal**2*voxelsize_axial,0)))
    for i,row in enumerate(drusenList[::-1]): 
        metricsTree.insert("", END,text=lastIdx-i, values=row[1])
    metricsTree.pack()

    return drusenmap