""" 
Evaluation: Module to calculate MSE and STDDEV of two volumes
------------------------------------------------------------------------
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

def calc_MSE_STDDEV( vol_test, vol_gt, seg_val):
    """ 
        Module to calculate MSE and STDDEV of two volumes
        
        Parameters
        ----------
        vol_test:ndarray
            test volume
        vol_gt: ndarray
            ground truth volume
        seg_val: scalar
            value of layer to segment (e.g., 255= Bruchs Membrane)
            
        Returns
        ---------
        MSE: scalar, float 
            mean standard error over volume
        SDEV: scalar float
            average standard deviation over volume
    """
    gt = []
    test = []
    for z in range(vol_test.shape[0]):
        for x in range (vol_test.shape[2]):
            try:
                gt_val = np.where(vol_gt[z][:,x] == seg_val)[0][-1]
                test_val = np.where(vol_test[z][:,x] == seg_val)[0][-1]
                
                gt.append(gt_val)
                test.append(test_val)
            except:
                continue

    diff_array = np.abs(np.asarray(gt) - np.asarray(test))
    return np.mean(diff_array),np.std(diff_array)