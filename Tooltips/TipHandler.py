""" 
TipHandler: Module to create tool tips
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
from tkinter import *

class CreateToolTip(object):
    """
        Creating a tool tip for a given TK Widget
        
        Parameters
        ----------
        object: tk widget
            the desired tk widget
    """
    def __init__(self, tkobject, text=""):
        """
            initializer
            
            set up handlers to enter and close label
            
            Parameters
            ----------
            self: object
                class
            tkobject: tkobject
                the tkobject
            text: string
                the text to set
        """
        self.tkobject = tkobject
        self.text = text
        self.id = None
        self.waittime = 2000
        self.tkobject.bind("<Enter>", self.enter)
        self.tkobject.bind("<Leave>", self.close)
        
        
    def enter(self, event=None):
        """
            Create tool tip on focus.
            
            Parameters
            ----------
            self: object
                class
            event: event
                event          
        """
        x = y = 0
        x, y, _, _ = self.tkobject.bbox("insert")
        
        x += self.tkobject.winfo_rootx() + 25
        y += self.tkobject.winfo_rooty() + 30

        self.level = Toplevel(self.tkobject)

        self.level.wm_overrideredirect(True)
        self.level.wm_geometry("+%d+%d" % (x, y))
        ttip_label = Label(self.level, text=self.text, justify='left', background='#A0A0A0',foreground='#000000', relief='solid', borderwidth=2, font=("Times 11"))
        ttip_label.pack(ipadx=1)
        
    def close(self, event=None):
        if self.level:
            self.level.destroy()
            
def createToolTipDictionary():
    """
        Create a dictionary for setting tool tips.
        
        This could be just rewritten for translation different languages.
    """
    dict = {}
    #Data Menu
    dict['loadRaw'] = "Load OCT volume and flatten it in slow-scan direction.\n\nMake sure that the volume has not been flattened yet.\nFlattening is done for each B-scan as follows:\n1) Extracting the RPE\n2) Ransac fitting RPE\n3) Shifting center of RPE to B-scan center"
    dict['loadMerged'] ="Load OCT volume.\n\nThis function does not flatten the volume in slow-scan direction.\nMerged data processed by the motion compensation algorithms\nshould be loaded with this function (e.g., Kraus et al)."
    dict['loadOCTA'] ="Load OCT Angiography volume.\n\nThis function just loads a volume without any processing.\nIf special modifications for OCTA have to be done,\nModify this handler."
    dict['loadSeg'] ="Load complete segmentation.\n\nLoad a tif(f) file segmentation fitting the OCT volume."
    dict['resetSeg'] ="Reset the segmentation."
    dict['saveInpainted'] ="Save Segmentation inpainted in volume.\n\nThe three lines will be inpainted in the set colors\n(cf Segmentation menu) and saved as .tif file."
    dict['saveLines'] ="Save Segmentation lines in empty volume.\n\nThe three lines will be stored in an empty .tif volume."
    dict['saveHM'] ="Save the heatmap as image.\n\nStore the map currently shown on the top right as image."
    
    #image Enhancement Menu
    dict['contrast'] ="Set Brightness and Contrast for visualization.\n\nBrightness and contrast can be independently set to enhance the visual output."
    dict['thresholdToOne'] ="Set grayscale or binary output."
    dict['threshold'] ="Activate and reset global thresholding.\n\nConfigure a certain threshold by varying the slider.\nThe button above can be used to configure binary\nor grayscale output."
    dict['colormap'] ="Change the heatmap's colormap.\n\nIt is also possible to integrate own colormaps or existing ones.\nCheck module: matplotlib.cm for more information."
    
    #Segmentation Menu
    dict['manRefine'] ="Open manual refinement frame.\n\nThe correction is based on the algorithm of Stromer et al. but can be individually adapted."
    dict['threeLayer'] ="Run automatic three layer segmentation (ILM, RPE, BM).\n\nThe segmentation algorithm and can be described as follows:\n1) Flattening to RPE\n2) Calculating graph weights\n3) Shortest path search\n4) Propagation for neighboring B-scans.\n\nMore details in Stromer et al.'s article at BOE 2019."
    dict['RPERefine'] ="Run RPE refinement algorithm.\n\nNote: This step needs a prior segmentation of Bruch's Membrane!\nA graph-cut is run similar to the three layer algorithm.\nAdaptions:\ni) Extended neighborhood for higher accuracy (9 neighbors)\nii) Only area above Burch's Membrane allowed to exclude CC structure"
    dict['segInpaint'] ="Inpaint segmentation lines into volume.\n\nThe three segmentation lines will be inpainted with the set colors."
    dict['filInpaint'] ="Inpaint filtered segmentation lines into volume.\n\nBased on the thresholded heatmap, the sensitivity through thresholding\nwill be propagated to create new segmentation lines for RPE and BM."
    dict['resetInpaint'] ="Reset inpainted segmentation.\n\nDelete segmentation lines from visualization."
    dict['colorILM'] ="Color for inpainting Inner Limiting Membrane (ILM).\n\nFormat: Hexadecimal RGB = #RRGGBB (R=red, G=green, B=blue).\nPress Return to propagate."
    dict['colorRPE'] ="Color for inpainting Retinal Pigment Epithelium (RPE).\n\nFormat: Hexadecimal RGB = #RRGGBB (R=red, G=green, B=blue).\nPress Return to propagate"
    dict['colorBM'] ="Color for inpainting Bruch's Membrane (BM).\n\nFormat: Hexadecimal RGB = #RRGGBB (R=red, G=green, B=blue).\nPress Return to propagate"
    dict['flattBM'] ="Flatten B-scans to Bruch's Membrane\n\nThe columns will be shifted such that the Bruch's Membrane\nlies on a straight line in the center.\nNote:This step needs a prior segmentation of Bruchs' Membrane!\n"
    dict['flattRPE'] ="Flatten B-scans to RPE.\n\nThe columns will be shifted such that the RPE\nlies on a straight line in the center.\nA curve is fitted through the RPE with a polynomial of a certain degree.\n"
    
    #Visualization menu
    dict['rpedcthick'] ="Generate and show RPEDC thickness map.\n\nThe distance between RPE and BM is calculated and the offset subtracted from the result."
    dict['severeness'] ="Based on the RPEDC thickness map, drusen are sorted and color-coded by their volume.\n\nThe pop-up windows shows the drusen's volume in descending order\n\nThe sensitivity can be varied by thresholding the heatmap."
    dict['subrpe'] ="Generate and show the sub-RPE slab.\n\nThe thickness and projection mode can be varied."
    dict['subbm'] ="Generate and show the sub-BM slab.\n\nThe thickness and projection mode can be varied."
    dict['ilmrpe'] ="Generate and show the ILM to RPE projection.\n\nThe projection mode can be varied."
    dict['slabthickness'] ="Select the thickness of the slab given in pixel.\nPress return to propagate."
    dict['projectionmode'] ="Select the projection mode:\nMean: average of axial direction\nMedian: median of axial direction\nMinmum: minimum of axial direction\nMaximum: maximum of axial direction"
    dict['vesselness'] ="Calculate a vessel probability map.\n\nThe probabilities are calculated by applying Frangi et al.'s Vesselness filter."
    dict['delineateGA'] ="Semi-automatically delineate GA regions.\n\nBy using a Grab-Cut algorithm, GA (or bright) lesions can be semi-automatically delineated.\nDraw a rectangle around a region-of-interest."
    dict['blender'] ="Inpaint the heatmap into the volume with a given percentage configured by the slider."
  
    #Tools menu
    dict['transversal'] ="Transversal resolution in microns."
    dict['axial'] ="Axial resolution in microns."
    dict['activateRes'] ="Activate resolution changes."
    #planes-canvas-labels
    dict['heatmapThresh'] ="Configure threshold value.\n\nUsed to set sensitivity of heatmap."
    dict['heatmapIntensity'] ="Show intensity.\n\nRight click on heatmap provides the value at the given pixel."
    dict['enfaceIntensity'] ="Show intensity.\n\nRight click on specific plane provides the value at the given pixel."
    
    #manual correction frame
    dict['manExplore'] ="Change back to volume explore mode.\n\nLeaving the manual correction frame. Segmentation will be transferred."
    dict['manLoadGT'] ="Load ground truth data.\n\nSelect ground truth data to show deviations to automatically segmented result."
    dict['manAutoEval'] ="Run automatic algorithm evaluation.\n\nDefine r,l, and spacing in auto_evaluation.txt. The evaluation then evaluates MSE and STDDEV\nfor the algorithm of the current volume. The algorithm uses the ground truth as input\nfor the corrected lines and sets an equal spacing.Example:\n1) Load Volume\n2) Run automatic segmentation\n3) Switch to Manual Refinement\n4) Load ground truth\n\nCheck publication of Stromer et al. for more information."
    dict['manResetSlice'] ="Reset manual correction of current slice and insert original segmentation."
    dict['manSaveRes'] ="Save the current segmentation."
    dict['manRunRefine'] ="Run manual refinement algorithm based on Stomer et al.'s algorithm.\n\n1) Correct a (set of) B-scan(s)\n2) Limit the search region by drawing a rectangle in the en face plane.\n3) Press this button to run the algorithm."
    dict['manCorrLayer'] ="Select Layer to be corrected.\n\nILM: Inner Limit Membrane\nRPE: Retinal Pigment Epithelium\nBM: Bruch's Membrane"
    return dict