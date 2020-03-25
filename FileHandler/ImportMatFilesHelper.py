'''
Created on 19.03.2020

@author: Daniel
'''
import scipy.io
import numpy as np

def loadmatfiles(filename):
    dict = scipy.io.loadmat(filename)
    vol = dict['images'].astype('int32')*255
    vol = np.swapaxes(np.asarray(vol), axis1 = 0, axis2 = 2)
    vol = np.swapaxes(np.asarray(vol), axis1 = 1, axis2 = 2)
    return vol

def loadmatfiles_segmentation(filename,layer_values):
    dict = scipy.io.loadmat(filename)
    vol = dict['images']
    vol_layers = dict['layerMaps']
    
    try: 
        seg_vol = np.zeros((vol.shape[2],vol.shape[0],vol.shape[1])).astype('float32')
    except Exception as e:
        print(e)
    
    vol_ilm = vol_layers[:,:,0].astype('float32')
    vol_rpe = vol_layers[:,:,1].astype('float32')
    vol_bm = vol_layers[:,:,2].astype('float32')
    
    for z in range (seg_vol.shape[0]):
        for x in range (seg_vol.shape[2]):
            if np.isfinite(vol_ilm[z,x]):
                seg_vol[z,vol_ilm[z,x].astype('int32'),x] = layer_values[0]
            if np.isfinite(vol_rpe[z,x]):
                seg_vol[z,vol_rpe[z,x].astype('int32'),x] = layer_values[1]
            if np.isfinite(vol_bm[z,x]):
                seg_vol[z,vol_bm[z,x].astype('int32'),x] = layer_values[2]
            
    return seg_vol.astype('int32')