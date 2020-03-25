""" 
ImportHandler: Module to import files
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
import os
from Visualization.Heatmap import calculateHeatmap
from Algorithms.SlowScanRegistration import RPEBasedRegistration
import scipy.io
from FileHandler import ImportMatFilesHelper

def loadVolume(directory, path=None, OCTA = False, merged = True, flattening_polynomial = 4):
    """ 
        Load volume from '.tif' stack
        
        Parameters
        ----------
        directory : string
            current path
        path : string, optional
            given path for skipping manual selectionc
        OCTA : bool, optional
            True if loading OCTA scan
        merged: bool
            volume merged
        flattening_polynomial: scalar
            degree for flattening
            
        Returns
        ----------
        initialdir : string
            new current path
        shifts : list of ndarrays
            the A-scan shift values for flattening: format [B-scan number][Column,y-coordinate]
        volume_original : ndarray
            unflattened volume
        vol_rgb: ndarray
            volume in rgb to visualize
    """    
    if(path is None):
        filename =  filedialog.askopenfilename(initialdir = directory,title = "Select file",filetypes = (("all files","*.*"),("tiff files","*.tiff"),("mat files","*.mat")))
    else:
        filename = path
    
    initialdir = os.path.split(filename)[0]
    f = open("options.dat", "w")
    f.write("dir="+initialdir)
    f.close()
    
    isMatFile = False
    if(os.path.splitext(filename)[1] == '.mat'):
        vol = ImportMatFilesHelper.loadmatfiles(filename).astype('float32')
        isMatFile = True
    else:
        vol = io.imread(filename).astype('float32')
    shifts = np.zeros((vol.shape[0],vol.shape[2]))

    try:
        if merged == False:
            vol, shifts = RPEBasedRegistration.runSlowScanRegistration(vol, 2)
    except Exception as e:
        print("Slow scan registration did not work")
    volume_original = vol.copy()
    
    vol = np.swapaxes(np.asarray(vol), axis1 = 0, axis2 = 1)
    vol_rgb = np.zeros((vol.shape[0],vol.shape[1],vol.shape[2],3)).astype('uint8')
    vol_b = (255*vol/65535.0).astype('uint8')
    vol_rgb[:,:,:,0] = vol_b[:,:,:]
    vol_rgb[:,:,:,1] = vol_b[:,:,:]
    vol_rgb[:,:,:,2] = vol_b[:,:,:]
    
    
    return initialdir, np.asarray(shifts), np.asarray(volume_original), vol_rgb, (isMatFile, filename)

def loadGroundtruthVolume(directory, shiftedValues):
    """ Load volume from '.tif' stack
    
    Parameters
    ----------
    directory : string
        current path
    shiftedValues : list of ndarrays
        the A-scan shift values for flattening: format [B-scan number][Column,y-coordinate]

    Returns
    ----------
    vol : string
        the groundtruth volume

    """    
    
    filename =  filedialog.askopenfilename(initialdir = directory,title = "Select file",filetypes = (("tif files","*.tif"),("all files","*.*")))
    
    vol = io.imread(filename).astype('int32')
    for z in range(vol.shape[0]):
        for x in range (vol.shape[2]):
            vol[z,:,x] = np.roll(vol[z,:,x], shiftedValues[z][x], axis= 0)
    
    return vol

def getOriginalRGBVolume(volume):
    """ Get original rgb volume
    
    Parameters
    ----------
    volume : ndarray
        volume
        
    Returns
    ----------
    vol_rgb : ndarray
        original rgb volume

    """ 
    vol_b = np.swapaxes(np.asarray(volume), axis1 = 0, axis2 = 1)
    vol_rgb = np.zeros((vol_b.shape[0],vol_b.shape[1],vol_b.shape[2],3)).astype('uint8')
    vol_b = (255*vol_b/65535.0).astype('uint8')
    
    vol_rgb[:,:,:,0] = vol_b[:,:,:]
    vol_rgb[:,:,:,1] = vol_b[:,:,:]
    vol_rgb[:,:,:,2] = vol_b[:,:,:]
    
    return vol_rgb

def loadCompleteSegmentation(directory, shifts, path=None):
    """ Load Segmentation
    
    Parameters
    ----------
    directory : string
        current path
    shifts : list of ndarrays
        the A-scan shift values for flattening: format [B-scan number][Column,y-coordinate]
        
    Optional
    ----------
    path : string, optional
        given path for skipping manual selection
        
    Returns
    ----------
    segmentation : ndarray
        the segmentation
    heatmap: ndarray
        calculated rpe/bruchs heatmap
    """ 
    if path is None:    
        filename_segmentation =  filedialog.askopenfilename(initialdir = directory,title = "Select file",filetypes = (("tif files","*.tif"),("all files","*.*")))
        segmentation = io.imread(filename_segmentation)
        f = open("options.dat", "w")
        f.write("dir="+os.path.split(filename_segmentation)[0])
        f.close()
    else:
        segmentation = io.imread(path)
    
    try:
        shifts = shifts.astype('int32')
        for z in range(segmentation.shape[0]):
            for x in range (segmentation.shape[2]):
                segmentation[z,:,x] = np.roll(segmentation[z,:,x],shifts[z][x], axis= 0)
    except Exception as e:
        print('No volume loaded',e)

    heatmap = calculateHeatmap(segmentation)
    return segmentation, heatmap

def loadCompleteSegmentationFromMat(isMatFile, shifts, path=None, layer_values=None):
    """ Load Segmentation
    
    Parameters
    ----------
    directory : string
        current path
    shifts : list of ndarrays
        the A-scan shift values for flattening: format [B-scan number][Column,y-coordinate]
        
    Optional
    ----------
    path : string, optional
        given path for skipping manual selection
        
    Returns
    ----------
    segmentation : ndarray
        the segmentation
    heatmap: ndarray
        calculated rpe/bruchs heatmap
    """ 
 
    segmentation = ImportMatFilesHelper.loadmatfiles_segmentation(isMatFile[1],layer_values).astype('float32')

    try:
        shifts = shifts.astype('int32')
        for z in range(segmentation.shape[0]):
            for x in range (segmentation.shape[2]):
                segmentation[z,:,x] = np.roll(segmentation[z,:,x],shifts[z][x], axis= 0)
    except Exception as e:
        print('No volume loaded',e)

    heatmap = calculateHeatmap(segmentation)
    return segmentation, heatmap