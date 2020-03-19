""" 
ExportHandler: Module to Export Files
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
from tkinter import filedialog
import numpy as np
from skimage import io
from Algorithms.Flattening.VolumeFlattening import unFlatten
import os 
from matplotlib import cm

def saveInpainted(directory, volumeXZ, shiftedValues):
    """ 
        Export lines inpainted in original volume as tif/tiff
    
        Parameters
        ----------
        directory : string
            current path
        volumeXZ : ndarray
            Data stack
        shiftedValues : scalar
            Flattening values for every A-scan
    """
    filename =filedialog.asksaveasfilename(initialdir = directory,title = "Save as ...",filetypes = (("all files","*.*"),("tiff","*.tiff"),("tif","*.tif")))
    f = open("options.dat", "w")
    f.write("dir="+os.path.split(filename)[0])
    f.close()
    
    if filename is None: 
        return volumeXZ

    if np.logical_not(filename[-4:] == '.tif' or filename[-5:] =='.tiff'):
        filename += '.tif'

    saveVolume = volumeXZ.copy()
    try:
        saveVolume = unFlatten(saveVolume , shiftedValues.astype('int32')).astype('uint8')
    except Exception as e:
        print('Volume could be empty:\n',e)
    io.imsave(filename, saveVolume)
    
def saveLineSegmentation(directory, segmentation, shiftedValues):
    """ 
        Export lines inpainted in empty volume as tif/tiff
    
        Parameters
        ----------
        directory : string
            current path
        segmentation : ndarray
            Data stack
        shiftedValues : scalar
            Flattening values for every A-scan
    """
    filename =filedialog.asksaveasfilename(initialdir = directory,title = "Save as ...",filetypes = (("all files","*.*"),("tiff","*.tiff"),("tif","*.tif")))
    if filename is None: 
        return

    if np.logical_not(filename[-4:] == '.tif' or filename[-5:] =='.tiff'):
        filename += '.tif'
        
    f = open("options.dat", "w")
    f.write("dir="+os.path.split(filename)[0])
    f.close()
    
    saveVolume = segmentation.copy()
    try:
        saveVolume = unFlatten(saveVolume, shiftedValues.astype('int32')).astype('uint8')
    except:
        print('Volume could be empty')
    io.imsave(filename, saveVolume)
    
def saveHeatmap(directory, heatmapForSaving, COLORMAP):
    """ 
        Export heatmap
        
        Parameters
        ----------
        directory : string
            current path
        heatmapForSaving : ndarray
            heatmap
        COLORMAP: string
            colormap settings 

    """
    cmap_std= cm.get_cmap(COLORMAP, np.max(heatmapForSaving))
    filename = filedialog.asksaveasfilename(initialdir = directory,title = "Save as ...",filetypes = (("all files","*.*"),("tiff","*.tiff"),("tif","*.tif"), ("png","*.png"),("jpeg","*.jpeg")))
    if filename is None: 
        return
    
    f = open("options.dat", "w")
    f.write("dir="+os.path.split(filename)[0])
    f.close()
    
    if '.tif' not in filename[-8:] and '.tiff' not in filename[-8:] and '.png' not in filename[-8:] and '.jpeg' not in filename[-8:]:
        filename += '.png'
    io.imsave(filename, cmap_std(heatmapForSaving))