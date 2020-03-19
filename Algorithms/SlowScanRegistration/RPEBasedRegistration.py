""" 
RPEBasedRegistration: Module to register the volume in slow scan direction
-----------------------------------------------------------------------------------
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
#sigma for Gaussian filter
sigma = (6.0,3.0)

def runSlowScanRegistration(volume, flattening_polynomial):
    """
        Flattening an OCT volume in slow-scan direction.
        
        Each B-scan is Gaussian smoothed and the index of the maximum 
        of each A-scan stored. The resulting pints are ransac-filtered
        and fitted by a polynomial of 2nd order. Finally the B-scan center
        index is extracted the entire B-scan shifted, such that the center 
        A-scan is located in the axial center.
        
    """
    #return volumes
    shiftedVolume = []
    shiftedValues = []
    
    im_center_y = int(volume[0].shape[0]//2)

    counter = 0
    prev_y = 0

    x_idx = np.arange(volume[0].shape[1]).astype(np.int32)
    #fill missing columns with neighbors
    volume = fillMissingColumns(volume).astype('float32')

    #iterate through slices
    for i,_slice in enumerate(volume):
        #Gaussian blurring
        blurred = cv2.GaussianBlur(_slice, (0,0), sigmaX = sigma[0], sigmaY=sigma[1])
        # max of each A-scan
        y_idx = np.argmax(blurred, axis = 0)
        #outlier removal
        try:
            degree = flattening_polynomial
            x_idx, y_idx = ransac_fit(x_idx, y_idx, degree)
        except:
            try:
                degree = flattening_polynomial
                x_idx, y_idx = ransac_fit(x_idx, y_idx, degree) 
            except:
                y_idx = prev_y
        #get center y and translate volume
        fit_center_y = y_idx[y_idx.size//2]
        roll_value = np.zeros((y_idx.size)).astype('int32') + (im_center_y-fit_center_y)
        new_slice = np.roll(_slice,roll_value[0], axis= 0)
        #store B-scan and shifts
        shiftedValues.append(roll_value)
        shiftedVolume.append(new_slice)
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