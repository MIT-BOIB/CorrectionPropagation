""" 
VolumeFlattening: Module to flatten and un-flatten OCT data based on RPE
--------------------------------------------------------------------------------
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
import cv2
from sklearn import linear_model
from sklearn.preprocessing import PolynomialFeatures
from FileHandler import ImportHandler
#sigma for Gaussian filter
sigma = (4.0,2.0)

def runFlattening(volume, flattening_polynomial):
    """
        Flattening an OCT volume based on RPE. 
        
        Compare S. Chiu's paper.
        
        Algorithm:
        For every volume slice...
        1) Gaussian Smoothing
        2) Searching brightest pixel in every column, and only allow +/-3 pixels distance from prior slice
        3) Ransac and curve fit of degree 4
        4) Shift columns and store shifts
        
        Parameters
        ----------
        volume: numpy array 2D/3D array
            input oct volume_rgb
        flattening_polynomial = scalar
            input for polynomial fitting - order of the curve
        
        Return
        ------
        shiftedVolume: list of 2D numpy arrays
            flattened volume_rgb slices
        shiftedValues: list of list
            the shifted values of each column - if the y-th column of slice z was shifted/translated by t,
            the translation distance t can be accessed by shiftedValues[z][y].
        
    """
    #return volumes
    shiftedVolume = []
    shiftedValues = []
    im_center_y = int(volume[0].shape[0]//2)

    counter = 0
    prev_y = 0
    
    x_idx = np.arange(volume[0].shape[1]).astype(np.int32)
    #fill missing columns with neighbors
    volume = fillMissingColumns(volume)
    #iterate through slices
    for i,_slice in enumerate(volume):
        #blurring - kernelwidth computed from sigmas (0,0)
        blurred = cv2.GaussianBlur(_slice, (0,0), sigmaX = sigma[0], sigmaY=sigma[1])
        y_idx = 0
        # get max of colums
        if(i == 0):
            y_idx = np.argmax(blurred, axis = 0)
        else:
            for x in range (blurred.shape[1]):
                blurred[0:int(prev_y[x])-3,x] = 0
                blurred[int(prev_y[x])+3:blurred.shape[1],x] = 0
            y_idx = np.argmax(blurred, axis = 0)
        #ransac outlier removal and polynomial fitting
        try:
            degree = flattening_polynomial
            x_idx, y_idx = ransac_fit(x_idx, y_idx, degree)
        except:
            try:
                degree = flattening_polynomial
                x_idx, y_idx = ransac_fit(x_idx, y_idx, degree) 
            except:
                y_idx = prev_y
            
        prev_y = y_idx
        #shifting volume
        shifts = [-int((y_idx_item - im_center_y)) for y_idx_item in y_idx]
        slicebuffer = np.zeros((_slice.shape)).astype('float32')
        for x in range (len(x_idx)):
            slicebuffer[:,x] = np.roll(volume[i][:,x],shifts[x], axis= 0)
        
        shiftedValues.append(shifts)
        shiftedVolume.append(slicebuffer)
        counter += 1

    return shiftedVolume, shiftedValues

def ransac_fit(x, y, degree_in):
    """
        Ransac fitting helper
        
        Parameters
        ----------
        x: ndarray
            data points x
        y: ndarray
            data points y
        degree_in: scalar
            degree of ransac
            
        Return
        ------
        x_i:ndarray
            data points x_i
        y_i: ndarray
            data points y_i
    """
    x_ = x.reshape((-1, 1))
    y_ = y.reshape((-1, 1))

    xi = np.linspace(min(x), max(x), len(x)).reshape((-1, 1))    

    m = linear_model.RANSACRegressor(linear_model.LinearRegression(), min_samples = 200, max_trials=100, random_state = 42)
    poly_2 = PolynomialFeatures(degree= degree_in,)
    x_2 = poly_2.fit_transform(x_)
    xi_2 = poly_2.fit_transform(xi)
    m.fit(x_2, y_)
    yi = m.predict(xi_2)
    
    return xi.astype(np.int32),yi.astype(np.int32)

def fillMissingColumns(volume):
    """
        This method fills up empty columns by its nearest neighbors.
        If columns are zero, RANSAC may fail on several slices. Furthermore, weights are not defined.
        
        Parameters
        ----------
        volume: list of ndarrays
            Input oct volume_rgb
            
        Return
        ------
        volume: list of ndarrays
            Column filled volume_rgb
    """

    #get all zero cols
    missing_cols = [(i,np.where(~_slice.any(axis=0))[0]) for i,_slice in enumerate(volume)]
    
    center = len(missing_cols)//2
    i = 1
    while True:
        if center - i >= 0:
            empty_cols = missing_cols[center - i][1]
            if empty_cols.size != 0:
                volume[center - i][:, empty_cols] = volume[center - i + 1][:, empty_cols]
            
        if center + i <len(missing_cols):
            empty_cols =  missing_cols[center + i][1]
            if empty_cols.size != 0:
                volume[center + i][:,empty_cols] = volume[center + i - 1][:,empty_cols]
        i += 1
        
        if center - i < 0 and center + i >= len(missing_cols):
            break
    
    return volume

def unFlatten(volume, shiftedValues):
    """
        Method to unflatten the slices
        
        Parameters
        ----------
        volume: ndarray
            flattened oct volume
        shiftedValues: list of lists
            the shifted values of each column - if the y-th column of slice z was shifted/translated by t,
            the translation distance t can be accessed by shiftedValues[z][y].
            
        Return
        ------
        volume: ndarray
            unflattened volume
    """
    #process all volume_rgb slices
    for i,_slice in enumerate(volume):
        for x in range (len(shiftedValues[i])):
            _slice[:,x] = np.roll(_slice[:,x],-shiftedValues[i][x], axis= 0)
    
    return volume

def applyFlattening(self, event=None):
    """
        Method called by GUI to apply RPE flattening and update planes
        
        Optional
        ----------
        event: event
            event

    """
    #reset BM flattening if active
    if self.flattenBM is False:
        applyFlatteningBM()
    
    if self.flatten is True:
        self.buttonFlattening.configure(text="Reset RPE Flattening", bg=self.btn_common_bg)
        #set polynomial
        try:
            self.flattening_polynomial = int(self.polynomial_var.get())
        except:
            self.polynomial_var.set(4)
            self.flattening_polynomial = 4
        #run algorithm
        vol, self.shiftedValues = runFlattening(self.volume_original, self.flattening_polynomial)
        self.shiftedValues = np.asarray(self.shiftedValues).astype('int32')
        # load as rgb
        self.volume = ImportHandler.getOriginalRGBVolume(np.asarray(vol))
        
        #-shiftedvalues flattens the image (note the minus)
        self.segmentation = unFlatten(self.segmentation,-self.shiftedValues)
        #update planes
        self.updateVolumeSliceXY(None)
        self.updateVolumeXZ()
        self.updateVolumeYZ()
        self.updateSegmentation()    
        self.flatten = False
    else:
        #reset flattening
        self.buttonFlattening.configure(text="Apply Flattening to RPE", bg=self.btnbackground)
        
        vol = self.volume_original.copy()
        self.volume = ImportHandler.getOriginalRGBVolume(np.asarray(vol))
        
        self.segmentation = unFlatten(self.segmentation, self.shiftedValues)
        
        self.updateVolumeSliceXY(None)
        self.updateVolumeXZ()
        self.updateVolumeYZ()
        self.updateSegmentation()
        self.flatten = True    
        self.shiftedValues = np.zeros((self.shiftedValues.shape[0],self.shiftedValues.shape[1])).astype('int32')
            
def applyFlatteningBM(self, event=None):
    """
        Apply Flattening to Bruch's Membrane Callback.
        
        Setting and Resetting Flattening to Bruch's Membrane.
        Updating planes and  Setting Buttons.
        
        Optional
        ----------
        event: event
            event
    """
    #reset RPE flattening if active
    if self.flatten is False:
        applyFlattening(self)
    
    if self.flattenBM is True:
        self.buttonFlatteningBM.configure(text="Reset Bruch's Flattening", bg=self.btn_common_bg)
        #run algorithm
        im_center_y = int(self.segmentation.shape[1]//2)  
        
        self.shiftedValues = np.zeros((self.shiftedValues.shape[0],self.shiftedValues.shape[1])).astype('int32')
        vol = self.volume_original.copy()
        for z in range(vol.shape[0]):
            coordinates_bruchs = np.where(self.segmentation[z].transpose() == self.BM_VALUE)[1]
            self.shiftedValues[z,:]= -coordinates_bruchs[:] + im_center_y
            for x in range (vol.shape[2]):
                vol[z,:,x] = np.roll(vol[z,:,x],self.shiftedValues[z][x], axis= 0)  
                
        self.volume = ImportHandler.getOriginalRGBVolume(np.asarray(vol)) 
        #-shiftedvalues flattens the image (note the minus)
        self.segmentation = unFlatten(self.segmentation, -self.shiftedValues)
        #update planes
        self.updateVolumeSliceXY(None)
        self.updateVolumeXZ()
        self.updateVolumeYZ()
        self.updateSegmentation()    
        self.flattenBM = False
    else:
        #reset flattening
        self.buttonFlatteningBM.configure(text="Flatten to Bruch's Membrane", bg=self.btnbackground)
        
        vol = self.volume_original.copy()
        self.volume = ImportHandler.getOriginalRGBVolume(np.asarray(vol))

        self.segmentation = unFlatten(self.segmentation, self.shiftedValues)
        
        self.updateVolumeSliceXY(None)
        self.updateVolumeXZ()
        self.updateVolumeYZ()
        self.updateSegmentation()
        self.flattenBM = True  
        self.shiftedValues = np.zeros((self.shiftedValues.shape[0],self.shiftedValues.shape[1])).astype('int32') 