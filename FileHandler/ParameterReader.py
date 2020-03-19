""" 
ParameterReader: Module to read Parameters from text into dictionaries
-------------------------------------------------------------------------
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
import os
import sys
import numpy as np

def readGlobalParameters(self):
    """
        Reading global parameters from parameters.txt
        
        This allows changes of parameters.
        
        Parameters
        ----------
        self: object
            Framework
    """
    #read options to receive stored last path of operation
    with open(os.path.dirname(os.path.abspath(sys.argv[0]))+"\\options.dat") as f:
        for line in f:
            if '#' in line:
                continue
            if 'dir=' in line: 
                self.initialdir = line.split('=')[1]
    rpe_found = False
    ilm_found = False
    bm_found = False
    self.SMALLWINDOW=False
    with open(os.path.dirname(os.path.abspath(sys.argv[0]))+"\\parameters.txt") as f:
        for line in f:
            if '#' in line:
                continue
            if 'RPE_VALUE' in line: 
                rpe_found = True
                self.RPE_VALUE=int(line.split('=')[1])
            elif 'ILM_VALUE' in line: 
                ilm_found = True
                self.ILM_VALUE=int(line.split('=')[1])
            elif 'BM_VALUE' in line: 
                bm_found = True
                self.BM_VALUE=int(line.split('=')[1])
            elif 'SMALLWINDOW' in line: 
                if 'True' in line or 'true' in line or 'TRUE' in line:
                    self.SMALLWINDOW=True
    if not bm_found:
        print('No BM Value found! Assuming 255!')
        self.BM_VALUE=255
    if not ilm_found:
        print('No ILM Value found!Assuming 127!')
        self.ILM_VALUE=127
    if not rpe_found:
        print('No RPE Value found!Assuming 64!')
        self.RPE_VALUE=64
        
def readThreeLayerDict():
    """
        Reading parameters for Three layer segmentation from parameters.txt
        
        This allows changes of parameters during runtime.
        
        Returns
        ----------
        dictParameters: dictionary
            Parameters dictionary
    """
    dictParameters={}
    rpe_found = False
    ilm_found = False
    bm_found = False
    with open(os.path.dirname(os.path.abspath(sys.argv[0]))+"\\parameters.txt") as f:
        for line in f:
            if '#' in line:
                continue
            if 'RPE_VALUE' in line: 
                rpe_found = True
                dictParameters['RPE_VALUE']=int(line.split('=')[1])
            elif 'ILM_VALUE' in line: 
                ilm_found = True
                dictParameters['ILM_VALUE']=int(line.split('=')[1])
            elif 'BM_VALUE' in line: 
                bm_found = True
                dictParameters['BM_VALUE']=int(line.split('=')[1])
            elif 'AUTO_ILM_BF_ENFACE' in line:
                dictParameters['AUTO_ILM_BF_ENFACE']=np.array(line.split('=')[1].split(',')).flatten().astype('int32')
            elif 'AUTO_ILM_BF_BSCAN' in line:
                dictParameters['AUTO_ILM_BF_BSCAN']=np.array(line.split('=')[1].split(',')).flatten().astype('int32')
            elif 'AUTO_RPE_BF_ENFACE' in line:
                dictParameters['AUTO_RPE_BF_ENFACE']=np.array(line.split('=')[1].split(',')).flatten().astype('int32')
            elif 'AUTO_RPE_BF_BSCAN' in line:
                dictParameters['AUTO_RPE_BF_BSCAN']=np.array(line.split('=')[1].split(',')).flatten().astype('int32')
    if not bm_found:
        print('No BM Value found! Assuming 255!')
        dictParameters['BM_VALUE'] = 255
    if not ilm_found:
        print('No ILM Value found!Assuming 127!')
        dictParameters['ILM_VALUE'] = 127
    if not rpe_found:
        print('No RPE Value found!Assuming 64!')
        dictParameters['RPE_VALUE'] = 64
    return dictParameters
 
def readRPERefinementDict():
    """
        Reading  parameters for RPE Refinement from parameters.txt
        
        This allows changes of parameters during runtime.
        
        Returns
        ----------
        dictParameters: dictionary
            Parameters dictionary
    """
    dictParameters = {}
    rpe_found = False
    ilm_found = False
    bm_found = False
    with open(os.path.dirname(os.path.abspath(sys.argv[0]))+"\\parameters.txt") as f:
        for line in f:
            if '#' in line:
                continue
            if 'RPE_VALUE' in line: 
                rpe_found = True
                dictParameters['RPE_VALUE']=int(line.split('=')[1])
            elif 'ILM_VALUE' in line: 
                ilm_found = True
                dictParameters['ILM_VALUE']=int(line.split('=')[1])
            elif 'BM_VALUE' in line: 
                bm_found = True
                dictParameters['BM_VALUE']=int(line.split('=')[1])
            elif 'REF_RPE_BF_ENFACE' in line:
                dictParameters['REF_RPE_BF_ENFACE']=np.array(line.split('=')[1].split(',')).flatten().astype('int32')
            elif 'REF_RPE_BF_BSCAN' in line:
                dictParameters['REF_RPE_BF_BSCAN']=np.array(line.split('=')[1].split(',')).flatten().astype('int32')
    if not bm_found:
        print('No BM Value found! Assuming 255!')
        dictParameters['BM_VALUE'] = 255
    if not ilm_found:
        print('No ILM Value found!Assuming 127!')
        dictParameters['ILM_VALUE'] = 127
    if not rpe_found:
        print('No RPE Value found!Assuming 64!')
        dictParameters['RPE_VALUE'] = 64
    return dictParameters
 
def readManRefParameters():
    """
        Reading parameters from parameters.txt into dictionary for manual refinement.
        
        This allows changes of parameters in runtime.
        
        Returns
        ----------
        dictParameters: dictionary
            Parameters dictionary
    """
    rpe_found = False
    ilm_found = False
    bm_found = False
    dictParameters = {}

    with open(os.path.dirname(os.path.abspath(sys.argv[0]))+"\\parameters.txt") as f:
        for line in f:
            if '#' in line:
                continue
            if 'RPE_VALUE' in line: 
                rpe_found = True
                dictParameters['RPE_VALUE']=int(line.split('=')[1])
            elif 'ILM_VALUE' in line: 
                ilm_found = True
                dictParameters['ILM_VALUE']=int(line.split('=')[1])
            elif 'BM_VALUE' in line: 
                bm_found = True
                dictParameters['BM_VALUE']=int(line.split('=')[1])
            elif 'MAN_ILM_BF_BSCAN' in line: 
                dictParameters['MAN_ILM_BF_BSCAN']=np.array(line.split('=')[1].split(',')).flatten().astype('int32')
            elif 'MAN_RPE_BF_BSCAN' in line: 
                dictParameters['MAN_RPE_BF_BSCAN']=np.array(line.split('=')[1].split(',')).flatten().astype('int32')
            elif 'MAN_ILM_MEDIAN' in line:  
                dictParameters['MAN_ILM_MEDIAN']=int(line.split('=')[1])
            elif 'MAN_RPE_MEDIAN' in line: 
                dictParameters['MAN_RPE_MEDIAN']=int(line.split('=')[1])
            elif 'MAN_BM_MEDIAN' in line: 
                dictParameters['MAN_BM_MEDIAN']=int(line.split('=')[1])
            elif 'MAN_ILM_THICKNESS' in line: 
                dictParameters['MAN_ILM_THICKNESS']=int(line.split('=')[1])
            elif 'MAN_RPE_THICKNESS' in line: 
                dictParameters['MAN_RPE_THICKNESS']=int(line.split('=')[1])
            elif 'MAN_BM_THICKNESS' in line: 
                dictParameters['MAN_BM_THICKNESS']=int(line.split('=')[1])
    if not bm_found:
        print('No BM Value found! Assuming 255!')
        dictParameters['BM_VALUE']  = 255
    if not ilm_found:
        print('No ILM Value found!Assuming 127!')
        dictParameters['ILM_VALUE'] = 127
    if not rpe_found:
        print('No RPE Value found!Assuming 64!')
        dictParameters['RPE_VALUE'] = 64
        
    return dictParameters