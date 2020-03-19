""" 
Framework: Core class to load initial GUI
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
import PIL.Image, PIL.ImageTk, PIL.ImageEnhance
from matplotlib import cm
from Visualization import Heatmap
from Visualization import Projections 
from Visualization import DelineateGA
from FileHandler import ExportHandler, ImportHandler
from Algorithms.AutomaticSegmentation import ThreeLayerSegmentation
from Algorithms.RefinementRPE import RPERefinementAlgorithm
from Algorithms.Flattening.VolumeFlattening import unFlatten
from GUI_Classes.OptionsPane import createOptionMenu
from GUI_Classes.CanvasXZ import createCanvasXZ
from GUI_Classes.CanvasXY import createCanvasXY
from GUI_Classes.CanvasYZ import createCanvasYZ
from GUI_Classes.CanvasHeatmap import createCanvasHeatmap
from GUI_Classes.FrameManualCorrection import createManualCorrectionFrame
from Tooltips.TipHandler import createToolTipDictionary
from FileHandler.ParameterReader import readGlobalParameters

class GUI:
    """ 
        Basic GUI Class to load PRLEC Framework
    """   
    #GUI colors
    btnbackground    = '#303030' #Button background
    btnforeground    = '#FFFFFF' #Button font
    btn_common_bg    = '#307130' #Button active
    canvasbackground = '#070707' #Canvas background
    canvasforeground = '#FAFAFA' #Canvas Font
    hlghtbg          = "#32E3E3" #Canvas highlighting
    crosshair_color  = "#32E3E3" #Crosshair color
    
    def __init__(self):
        """ 
            Initializing the GUI 
        """    
        #empty basic volumes for initial loader
        self.volume          = np.zeros((500,500,500)).astype('uint8') 
        self.volume_original = np.zeros((500,500,500)).astype('uint8')
        self.segmentation    = np.zeros((500,500,500)).astype('uint8')
        self.heatmap_rpe_bruchs  = np.zeros((500,500)).astype('uint8')
        self.heatmap_original    = np.zeros((500,500,3)).astype('uint8')
        self.metric_map          = np.zeros((500,500)).astype('uint8')
    
        self.ttip_dict = createToolTipDictionary()
        #heatmap minimum
        self.heatmap_min = 0
        '''
            Insert new i in parameters.txt
        '''
        readGlobalParameters(self)
        #GUI fonts
        if self.SMALLWINDOW is True:
            self.headerFont = ("Times 10 bold") #headlines
            self.textFont   = ("Times 9")      #text
            self.btnFont    = ("Times 9")      #Buttons
            self.statusFont = ("Times 9")      #statustext
            self.clrbarFont = ("Times 9")      #colorbar
        else:
            self.headerFont = ("Times 14 bold") #headlines
            self.textFont   = ("Times 13")      #text
            self.btnFont    = ("Times 13")      #Buttons
            self.statusFont = ("Times 12")      #statustext
            self.clrbarFont = ("Times 11")      #colorbar
        #Initial Heatmap
        self.heatmapType = 'Drusen'
        #Canvas Zoom variables
        self.scale_delta = 0.1
        self.scale_img = 1.0
        #Canvas Contrast variables
        self.contrast_value = 1.0
        self.brightness_value = 1.0
        #Flattening polynomial
        self.flattening_polynomial = 4
        #Rect for delineating GA
        self.rect = (0,0,0,0)
        #slow scan flattening
        self.shiftedValuesSlowScan = None
        #Variables to store buttons initial values
        self.blender = False
        self.contrast = False
        self.threshold_to_one = False
        self.threshold_activate = True
        self.flatten = True
        self.flattenBM = True
        self.projMode = 'mean'
        self.mover_init = True
        self.corr_active= False
        self.fromCorrMenu = False
        #blender job
        self._job_Blender = None
        #initial colormap
        self.COLORMAP = 'viridis'
        #Contours for GA
        self.contourlist = []
        #initial colors for segmentation inpainting
        self.colorRPE = StringVar()
        self.colorRPE.set("#FFBB00")
        self.colorBM = StringVar()
        self.colorBM.set("#BB0000")
        self.colorILM = StringVar()
        self.colorILM.set("#23CFCF")
        #resolutions
        self.tv_res = DoubleVar()
        self.tv_res.set(12)
        self.axial_res = DoubleVar()
        self.axial_res.set(4.5)
        #Visualization initially RPEDC
        self.varVis = IntVar()
        self.varVis.set(1)
    
    def start_explore_mode(self,master):
        """ 
            Start Explore Volume Mode
            
            Parameters
            ----------
            master: TK object
                TK object for setting up frame
        """
        #delete segmentation when coming from manual refinement
        if hasattr(self,'segmentation_original'):
            del self.segmentation_original
            
        #initial call
        if(not hasattr(self,'master')):
            self.master = master
            #geometry size
            if(self.SMALLWINDOW == True):  
                self.master.geometry('800x650')
            else:
                self.master.geometry('1200x850')
            #titlebar
            self.master.title("PRLE - OCT Processing Tool")
            self.master.state('zoomed')
            #binding arrow keys to zoom into planes
            self.master.bind("<Up>", self.keyPressedUp)
            self.master.bind("<Down>", self.keyPressedDown)
            self.master.configure(background=self.canvasbackground)
        #delete correction plane
        if(hasattr(self,'FrameSegCorrection')):   
            self.FrameSegCorrection.destroy()
        #Create menu's and planes for investigations
        createOptionMenu(self)
        self.fromCorrMenu = False
        createCanvasXY(self)
        createCanvasHeatmap(self)
        createCanvasXZ(self)
        createCanvasYZ(self)
        
    def exploremode(self):
        """  
            Starting Volume Explorer
        """
        #delete correction frame and set flattening
        if hasattr(self,"correctionFrame"):
            self.correctionFrame.destroy()
        if(hasattr(self,"shiftedValues") and np.count_nonzero(self.shiftedValues) != 0):
            self.flattenBM = False 
        self.corr_active= False
        self.start_explore_mode(self.master)
        
    def correctionmode(self):
        """  
            Starting Manual Refinement Mode
            Note: Volume will be flattened to BM for Correction
        """
        self.statusText.set("Switching to manual refinement...")
        #Check if volume is already flat
        if(self.flattenBM is False):
            self.start_correction_mode()
            self.corr_active= True
        else: 
            #Volume not flattened - flatten to BM
            if(np.count_nonzero(self.segmentation) != 0):
                im_center_y = int(self.segmentation.shape[1]//2)  
                self.shiftedValues = self.shiftedValues.astype('int32')
                vol = self.volume_original.copy()
                
                for z in range(vol.shape[0]):
                    coordinates_bruchs = np.where(self.segmentation[z].transpose() == self.BM_VALUE)[1]
                    self.shiftedValues[z,:]= -coordinates_bruchs[:] + im_center_y
                    for x in range (vol.shape[2]):
                        vol[z,:,x] = np.roll(vol[z,:,x],self.shiftedValues[z][x], axis= 0)  
                        
                self.volume = ImportHandler.getOriginalRGBVolume(np.asarray(vol))   
                #note the minus: -shiftedvalues flattens the image
                self.segmentation = unFlatten(self.segmentation,-self.shiftedValues)
                #update volumes
                self.updateVolumeSliceXY(None)
                self.updateVolumeXZ()
                self.updateVolumeYZ()
                self.updateSegmentation() 
                self.flattenBM = False
                
                #start correction mode
                self.start_correction_mode()
                self.corr_active= True
            #Volume not flattened and no base-segmentation exists
            else:
                self.segmentation = np.zeros((self.volume.shape[1],self.volume.shape[0],self.volume.shape[2])).astype('uint8')
                im_center_y = int(self.segmentation.shape[1]//2)  
                
                for z in range(self.segmentation.shape[0]):
                    self.segmentation[z,im_center_y,:] = self.BM_VALUE
                    
                self.shiftedValues = self.shiftedValues.astype('int32')
                vol = self.volume_original.copy()
                for z in range(vol.shape[0]):
                    coordinates_bruchs = np.where(self.segmentation[z].transpose() == self.BM_VALUE)[1]
                    self.shiftedValues[z,:]= -coordinates_bruchs[:] + im_center_y
                    for x in range (vol.shape[2]):
                        vol[z,:,x] = np.roll(vol[z,:,x],self.shiftedValues[z][x], axis= 0)  
                        
                self.volume = ImportHandler.getOriginalRGBVolume(np.asarray(vol))
                #note the minus: -shiftedvalues flattens the image
                self.segmentation = unFlatten(self.segmentation,-self.shiftedValues)
                #update volumes
                self.updateVolumeSliceXY(None)
                self.updateVolumeXZ()
                self.updateVolumeYZ()
                self.updateSegmentation() 
                self.flattenBM = False
                #start correction mode
                self.start_correction_mode()
                self.corr_active= True
        
        self.fromCorrMenu = True
    
    def saveInpainted(self):
        """  
            Saving inpainted volume as tif
            A window will be opened to select the path
        """
        self.statusText.set("Saving inpainted volume...")
        self.disableDataButtons()
        try:
            self.inpaintOriginalSegmentation()
            
            if self.shiftedValuesSlowScan is None:
                ExportHandler.saveInpainted(self.initialdir, self.volumeXZ, self.shiftedValues)
            else:
                ExportHandler.saveInpainted(self.initialdir, self.volumeXZ, self.shiftedValues+self.shiftedValuesSlowScan)
            
            self.statusText.set("Inpainted segmentation saved!")   
        
        except Exception as e:
            self.statusText.set("Attention! Segmentation not saved!")
            print("Error:",e)
        
        self.enableDataButtons()
    
    def saveLineSegmentation(self):
        """  
            Saving line segmentation as tif
            A window will be opened to select the path
        """
        self.statusText.set("Saving line segmentation...")
        self.disableDataButtons()
        
        try:
            
            if self.shiftedValuesSlowScan is None:
                ExportHandler.saveLineSegmentation(self.initialdir, self.segmentation, self.shiftedValues)    
            else:
                ExportHandler.saveLineSegmentation(self.initialdir, self.segmentation, self.shiftedValues+self.shiftedValuesSlowScan)    
                self.statusText.set("Line segmentation saved!")
        
        except Exception as e:
            self.statusText.set("Attention! Error while saving.")
            print('Error:',e)
        self.enableDataButtons()
        
    def saveHeatmap(self):
        """  
            Saving heatmap
            A window will be opened to select the path
        """
        self.statusText.set("Saving heatmap...")
        self.disableDataButtons()
       
        if self.heatmapType is 'minILMRPEProjection' or self.heatmapType is 'MarkGaByRect' or self.heatmapType is 'subRPESlab' or self.heatmapType is 'Vesselness':
            self.COLORMAP = 'gray'
        try:
            ExportHandler.saveHeatmap(self.initialdir, self.heatmapForSaving,self.COLORMAP)
            self.statusText.set("Heatmap saved!")
        except Exception as e:
            self.statusText.set("Attention: Heatmap not saved!")
            print('Error:',e)
        
        self.enableDataButtons()
    
    def loadVolumeOCTA(self):
        """   
            Load OCTA volume calling loadVolume
            A window will be opened to select the path
        """   
        try:
            self.loadVolume(path=None, OCTA = True)
            self.statusText.set("OCTA volume loaded!")
        except Exception as e:
            self.statusText.set("Error while loading volume!")
            print("Error:",e)
        
    def loadVolumeUnmerged(self):
        """   
            Load Raw OCT volume calling loadVolume with parameter merged = False.
            A window will be opened to select the path. The slow-scan direction
            will be registered to generate a smoother volume. For more information
            on the registration, please look at:
            Algorithms -> SlowScanRegistration -> runSlowScanRegistration
        """ 
        try:
            self.loadVolume(path=None, merged = False, OCTA = False)
            self.statusText.set("Volume loaded!")
        except Exception as e:
            self.statusText.set("Error while loading volume!")
            print("Error:",e)
    
    def disableDataButtons(self):
        """
            Disable Buttons to ensure that no multiple loading processes are triggered
        """
        self.btnLoadUnmerged.configure(state=DISABLED)
        self.btnLoadMerged.configure(state=DISABLED)
        self.btnLoadOCTA.configure(state=DISABLED)
        self.btnLoadSeg.configure(state=DISABLED)
        self.btnResetSeg.configure(state=DISABLED)
        self.btnSaveSegInpainted.configure(state=DISABLED)
        self.btnSaveSegLines.configure(state=DISABLED)
        self.btnSaveHM.configure(state=DISABLED)
        
    def enableDataButtons(self):
        """
            Enable Buttons after loading a volume
        """
        self.btnLoadUnmerged.configure(state=NORMAL)
        self.btnLoadMerged.configure(state=NORMAL)
        self.btnLoadOCTA.configure(state=NORMAL)
        self.btnLoadSeg.configure(state=NORMAL)
        self.btnResetSeg.configure(state=NORMAL)
        self.btnSaveSegInpainted.configure(state=NORMAL)
        self.btnSaveSegLines.configure(state=NORMAL)
        self.btnSaveHM.configure(state=NORMAL)
    
    def loadVolume(self, path=None, OCTA = False, merged = True):
        """ 
            Load volume calling ImportHandler.loadVolume and setting up planes
            A window will be opened to select the path
        """ 
        if merged is True:
            self.statusText.set("Loading merged OCT volume...")  
        elif OCTA is True:
            self.statusText.set("Loading OCTA volume...")  
        elif merged is False and OCTA is False:
            self.statusText.set("Loading raw volume...")  
        self.disableDataButtons()
        
        try:   

            if merged is False:
                self.initialdir,self.shiftedValuesSlowScan,self.volume_original,self.volume = ImportHandler.loadVolume(self.initialdir, path, OCTA, merged, self.flattening_polynomial)
                self.shiftedValues = np.zeros((self.shiftedValuesSlowScan.shape[0],self.shiftedValuesSlowScan.shape[1]))
            else:
                self.initialdir,self.shiftedValues,self.volume_original,self.volume = ImportHandler.loadVolume(self.initialdir, path, OCTA, merged, self.flattening_polynomial)
                self.shiftedValuesSlowScan = None
            
            self.resetSegmentation()
            #resetting the scale
            self.scale_img = 1.0
            #intitally volume is not flat
            self.flatten = True
        except Exception as e:
            self.enableDataButtons()
            self.statusText.set("Error while loading volume!")
            return 
        self.enableDataButtons()
        #updating scales for browsing through volume
        self.vol_sliderXY.destroy()
        if self.SMALLWINDOW is True:
            self.vol_sliderXY = Scale(self.FrameVolumeXY, from_=0, to=self.volume.shape[0]-1, length = self.widthXY//2, orient=HORIZONTAL, command=self.updateVolumeSliceXY, bg=self.canvasbackground, fg =self.canvasforeground)
        else:
            self.vol_sliderXY = Scale(self.FrameVolumeXY, from_=0, to=self.volume.shape[0]-1, length = self.widthXY, orient=HORIZONTAL, command=self.updateVolumeSliceXY, bg=self.canvasbackground, fg =self.canvasforeground)
        self.vol_sliderXY.grid(row=3,column=0)
        self.vol_sliderXY.set(self.volume.shape[0]//2)
        
        self.vol_sliderXZ.destroy()
        if self.SMALLWINDOW is True:
            self.vol_sliderXZ = Scale(self.FrameVolumeXZ, from_=0, to=self.volumeXZ.shape[0]-1, length = self.widthXZ//2,orient=HORIZONTAL, command=self.updateVolumeXZSlice, bg=self.canvasbackground, fg =self.canvasforeground)
        else:
            self.vol_sliderXZ = Scale(self.FrameVolumeXZ, from_=0, to=self.volumeXZ.shape[0]-1, length = self.widthXZ,orient=HORIZONTAL, command=self.updateVolumeXZSlice, bg=self.canvasbackground, fg =self.canvasforeground)
        self.vol_sliderXZ.grid(row=2,column=0)
        self.vol_sliderXZ.set(self.volumeXZ.shape[0]//2)
        
        self.vol_sliderYZ.destroy()
        if self.SMALLWINDOW is True:
            self.vol_sliderYZ = Scale(self.FrameVolumeYZ, from_=0, to=self.volumeYZ.shape[0]-1, length = self.widthYZ//2, orient=HORIZONTAL, command=self.updateVolumeYZSlice, bg=self.canvasbackground, fg =self.canvasforeground)
        else:
            self.vol_sliderYZ = Scale(self.FrameVolumeYZ, from_=0, to=self.volumeYZ.shape[0]-1, length = self.widthYZ, orient=HORIZONTAL, command=self.updateVolumeYZSlice, bg=self.canvasbackground, fg =self.canvasforeground)
        self.vol_sliderYZ.grid(row=2,column=0)
        self.vol_sliderYZ.set(self.volumeYZ.shape[0]//2)
        #updating volumes
        self.updateSegmentation()
        self.updateVolumeSliceXY(None)
        self.updateVolumeXZ()
        self.updateVolumeXZSlice(None)
        self.updateVolumeYZ()
        self.updateVolumeYZSlice(None)
        #loading segmentation - filename has to be complete_segmentation.tif
        try:
            self.loadCompleteSegmentation(self.initialdir+"//complete_segmentation.tif")
        except:
            print('Please load segmentation manually. No file found at:',self.initialdir+"//complete_segmentation.tif")
    
    def resetSegmentation(self):
        """ 
            Reset segmentation
        """ 
        self.statusText.set("Resetting segmentation...")  
        
        self.disableDataButtons()
        
        self.segmentation = np.zeros((self.volume.shape[1],self.volume.shape[0],self.volume.shape[2])).astype('uint8')
        self.heatmapType = 'Drusen'
        self.heatmap_rpe_bruchs = np.zeros((self.volume.shape[1],self.volume.shape[2])).astype('int32')
        self.updateSegmentation()
        self.enableDataButtons()
        self.statusText.set("Segementation reset done.")  
    
    def loadCompleteSegmentation(self, path=None):
        """ 
            Load segmentation volume calling ImportHandler.
        """ 
        self.statusText.set("Loading segmentation...")
        self.disableDataButtons()
        
        try:
            
            if self.shiftedValuesSlowScan is None:
                self.segmentation, self.heatmap_rpe_bruchs = ImportHandler.loadCompleteSegmentation(self.initialdir, self.shiftedValues, path)
            else:
                self.segmentation, self.heatmap_rpe_bruchs = ImportHandler.loadCompleteSegmentation(self.initialdir, self.shiftedValues+self.shiftedValuesSlowScan, path)   

            self.updateSegmentation()
            self.statusText.set("Segmentation loaded!")
        except Exception as e:
            self.statusText.set("Error while loading segmentation!")
            print("Error:", e)
        
        self.enableDataButtons()
 
    def updateSegmentation(self):
        """
            Updated Segmentation
            Different things happen for different heatmaps
        """
        if(self.heatmapType is 'Vesselness'):
            self.heat_slider_min.configure(from_=0, to=255)
            self.heat_slider_min.set(0)
            self.heatmap = self.vessel_original
            self.heatmap_max = 255
            self.updateHeatSlice(None)
            return
        elif(self.heatmapType is 'Drusen'):
            self.heatmap = self.heatmap_rpe_bruchs
            self.heatmapForSaving = self.heatmap
        elif(self.heatmapType is 'DrusenMetric'):
            self.heatmap = self.metric_map
            self.heatmapForSaving = self.heatmap
        elif(self.heatmapType is 'subRPESlab'):
            self.heatmapForSaving = self.heatmap.copy()
            
        self.heatmap_max = np.maximum(10,np.max(self.heatmap))
        self.heatmap_min = np.min(self.heatmap)
        
        self.label_top_HM.configure(text=' '+str(int(self.heatmap_max)))
        
        self.heat_slider_min.configure( to=self.heatmap_max, resolution=1)
        self.heat_slider_min.set(0)
        self.updateHeatSlice(None)
        
    def keyPressedUp(self, event):
        """
            Zooming into planes when key up pressed
        """
        self.scale_img = np.round(self.scale_img + self.scale_delta,1)
        self.updateVolumeSliceXY_Helper()
        self.updateVolumeYZSlice_Helper()
        self.updateVolumeXZSlice_Helper()
        self.updateHeatSlice_Helper()
    
    def keyPressedDown(self, event):
        """
            Zooming out of planes when key down pressed
        """  
        self.scale_img = np.round(self.scale_img - self.scale_delta,1)
        if(self.scale_img < 1.0):
            self.scale_img = 1.0

        self.updateVolumeSliceXY_Helper()
        self.updateVolumeYZSlice_Helper()
        self.updateVolumeXZSlice_Helper()
        self.updateHeatSlice_Helper()
        
    def getTopYMotion(self, *args):
        """
            Scrollbars vertical method
        """  
        self.canvasXY.yview(*args)
        self.canvas_heat.yview(*args)
        
    def getTopXMotion(self, *args):
        """
            Scrollbars horizontal method
        """    
        self.canvasXY.xview(*args)
        self.canvas_heat.xview(*args)
  
    def getBottomYMotion(self,*args):
        """
            Scrollbars vertical method
        """    
        self.canvasXZ.yview(*args)
        self.canvasYZ.yview(*args)
        
    def getBottomXMotion(self,*args):
        """
            Scrollbars horizontal method
        """      
        self.canvasXZ.xview(*args)
        self.canvasYZ.xview(*args)
           
    def updateVolumeSliceXY(self, event):
        """
            Updating XY plane image
            
            Parameters
            ----------
            event: event
                not used but necessary
        """    
        if self._job_volume_xy:
            self.master.after_cancel(self._job_volume_xy)
        self._job_volume_xy = self.master.after(50, self.updateVolumeSliceXY_Helper)
        
    def updateVolumeSliceXY_Helper(self):
        """
            Updating XY plane image helper
        """    
        self._job_volume_xy = None
        #setting up slice 
        try:
            self.slicexy = self.volume[self.vol_sliderXY.get()].astype('uint8')
        except:
            self.slicexy = self.volume[self.volume.shape[0]//2].astype('uint8')
            
        #setting up downscaled (self.SMALLWINDOW) or original size
        if(self.SMALLWINDOW == True):
            self.slice_pil_xy = PIL.Image.fromarray(self.slicexy[::2,::2])
        else:
            self.slice_pil_xy = PIL.Image.fromarray(self.slicexy)
        #zoom scale
        if self.scale_img is not 1.0:
            h,w = self.slice_pil_xy.size
            self.slice_pil_xy = self.slice_pil_xy.resize((int(h*self.scale_img), int(w*self.scale_img)), PIL.Image.ANTIALIAS)
        #contrast settings
        if self.contrast is True:
            self.slice_pil_xy = PIL.ImageEnhance.Contrast(self.slice_pil_xy).enhance(self.contrast_value)
            self.slice_pil_xy = PIL.ImageEnhance.Brightness(self.slice_pil_xy).enhance(self.brightness_value)
        #load image on canvas
        try:
            self.photoxy = PIL.ImageTk.PhotoImage(image = self.slice_pil_xy)
            self.canvasXY.itemconfig(self.image_on_canvas_xy, image = self.photoxy)
            self.canvasXY.configure(scrollregion = self.canvasXY.bbox("all"))
        except: 
            return
        #when blending active, blend
        self.blendSlices()

    def updateHeatSlice(self, event):
        """
            Updating Heatmap image
            
            Parameters
            ----------
            event: event
                not used but necessary
        """    
        if self._job_HM:
            self.master.after_cancel(self._job_HM)
        self._job_HM = self.master.after(50, self.updateHeatSlice_Helper)

    def updateHeatSlice_Helper(self):
        """
            Updating Heatmap image helper
        """  
        self._job_HM = None
        #load colormap
        self.cmap_own= cm.get_cmap(self.COLORMAP, self.heatmap_max)
        #Processing pipeline vesselness
        if(self.heatmapType is 'Vesselness'):
            #thresholding active
            if self.heat_slider_min.get() == 0.0:
                self.heatmapThreshed = self.vessel_original
            else: 
                self.heatmapThreshed = np.where(self.heatmap < self.heat_slider_min.get(), 0 , self.heatmap)
                self.heatmapThreshed = self.heatmapThreshed.astype('uint8')
            #configure colorbar
            cmap_std= cm.get_cmap('gray', 256)
            self.label_top_HM.configure(text=str(int(np.max(self.volume_original))))
            if(self.SMALLWINDOW == True):
                self.heatmapThreshedPil = PIL.Image.fromarray(self.heatmapThreshed[::2,::2])
                self.cbar_arr = np.array([range(0,self.heightHM//2)[::-1] for i in range(20)]).T.astype('uint8')
            else:
                self.heatmapThreshedPil = PIL.Image.fromarray(self.heatmapThreshed)    

                self.cbar_arr = np.zeros((self.heightHM, 20)).astype('uint8')       
                for i in range (256):
                    self.cbar_arr[i*2:i*2+2,:] = i
                self.cbar_arr = self.cbar_arr[::-1]   
        #Processing projections
        elif(self.heatmapType is 'minILMRPEProjection' or self.heatmapType is 'subRPESlab'):  
            #thresholding active
            if self.heat_slider_min.get() == 0.0:
                self.heatmapThreshed = self.heatmap.copy()
            else: 
                self.heatmapThreshed = np.where(self.heatmap < self.heat_slider_min.get(), 0 , self.heatmap)
                self.heatmapThreshed = self.heatmapThreshed.astype('uint8')
            #configure colorbar
            cmap_std= cm.get_cmap('gray', 256)
            self.label_top_HM.configure(text=str(int(np.max(self.volume_original))))
            if(self.SMALLWINDOW == True):
                self.heatmapThreshedPil = PIL.Image.fromarray(self.heatmapThreshed[::2,::2])
                self.cbar_arr = np.array([range(0,self.heightHM//2)[::-1] for i in range(20)]).T.astype('uint8')
            else:
                self.heatmapThreshedPil = PIL.Image.fromarray(self.heatmapThreshed)    

                self.cbar_arr = np.zeros((self.heightHM, 20)).astype('uint8')       
                for i in range (256):
                    self.cbar_arr[i*2:i*2+2,:] = i
                self.cbar_arr = self.cbar_arr[::-1]
        #Processing delineating GA
        elif(self.heatmapType is 'MarkGaByRect'):
            cmap_std= cm.get_cmap('gray', 256)
            self.label_top_HM.configure(text=str(int(np.max(self.volume_original))))
            if(self.SMALLWINDOW == True):
                self.heatmapThreshedPil = PIL.Image.fromarray(self.heatmapThreshed[::2,::2])
                self.cbar_arr = np.array([range(0,self.heightHM//2)[::-1] for i in range(20)]).T.astype('uint8')
            else:
                self.heatmapThreshedPil = PIL.Image.fromarray(self.heatmapThreshed)    

                self.cbar_arr = np.zeros((self.heightHM, 20)).astype('uint8')       
                for i in range (256):
                    self.cbar_arr[i*2:i*2+2,:] = i
                self.cbar_arr = self.cbar_arr[::-1]
        #Other mode
        else:
            #thresholding active
            try:
                self.heatmapThreshed = np.where(self.heatmap < self.heat_slider_min.get(), 0 ,self.heatmap).astype('uint8')
            except:
                self.heatmapThreshed = self.heatmap.astype('uint8')
            
            #configure colorbar
            cmap_std= cm.get_cmap(self.COLORMAP, 256)
            try:
                self.label_top_HM.configure(text=' '+str(np.max(self.heatmap)))
            except:
                return
            if(self.SMALLWINDOW == True):
                self.heatmapThreshedPil = PIL.Image.fromarray(self.cmap_own(self.heatmapThreshed[::2,::2], bytes=True))
                self.cbar_arr = np.array([range(0,self.heightHM//2)[::-1] for i in range(20)]).T.astype('uint8')
            else:
                self.heatmapThreshedPil = PIL.Image.fromarray(self.cmap_own(self.heatmapThreshed, bytes=True))      

                self.cbar_arr = np.zeros((self.heightHM, 20)).astype('uint8')       
                for i in range (256):
                    self.cbar_arr[i*2:i*2+2,:] = i
                self.cbar_arr = self.cbar_arr[::-1]
        #save original heatmap
        self.heatmap_original = self.heatmapThreshedPil.copy()
        
        #scale heatmap
        if self.scale_img is not 1.0:
            h,w = self.heatmapThreshedPil.size
            self.heatmapThreshedPil = self.heatmapThreshedPil.resize((int(h*self.scale_img), int(w*self.scale_img)), PIL.Image.ANTIALIAS)
            
        #update heatmap
        try:            
            self.photo_heat = PIL.ImageTk.PhotoImage(image = self.heatmapThreshedPil)
            self.photo_cbar = PIL.ImageTk.PhotoImage(image = PIL.Image.fromarray(cmap_std(self.cbar_arr, bytes=True)))    
            self.canvas_heat.itemconfig(self.image_on_canvas_heat, image = self.photo_heat)
            self.canvas_cbar_HM.itemconfig(self.image_on_canvas_cbar, image = self.photo_cbar)
        except: 
            return
        
        self.canvas_heat.configure(scrollregion = self.canvas_heat.bbox("all"))
        #update blender
        self.blendSlices()
 
    def blendSlices(self, event = None):
        """
            Updating Blender
            
            Parameters
            ----------
            event: event
                not used but necessary
        """    
        self._job_Blender = None
        #only valid blendings
        if(self.blender is False or self.heatmapType is 'Vesselness' or self.heatmapType is 'subRPESlab' or self.heatmapType is 'MarkGaByRect' or self.heatmapType is 'minILMRPEProjection'):
            return
        #update blender
        self.updateBlender(None)

    def updateBlender(self, event):
        """
            Updating Blender
            
            Parameters
            ----------
            event: event
                not used but necessary
        """    
        if self._job_Blender:
            self.master.after_cancel(self._job_Blender)
        self._job_XZ = self.master.after(50, self.updateBlenderHelper)
    
    def updateBlenderHelper(self):
        """
            Updating Blender Helper
            
            Parameters
            ----------
            event: event
                not used but necessary
        """    
        self._job_Blender = None
        #setting up image
        if(self.SMALLWINDOW == True):
            self.img_vol = self.slicexy[::2, ::2].astype('uint8')
        else:
            self.img_vol = self.slicexy.astype('uint8')    
        #contrast active
        if self.contrast is True:
            self.img_vol = PIL.ImageEnhance.Contrast(self.img_vol).enhance(self.contrast_value)
            self.img_vol = PIL.ImageEnhance.Brightness(self.img_vol).enhance(self.brightness_value)
        
        self.img_heatmap = np.asarray(self.heatmap_original)[:,:,0:3].astype('uint8')
        #add weighted images
        if(self.img_vol.size != self.img_heatmap.size):
            self.blended_img = np.zeros((self.img_heatmap.shape)).astype('uint8')
        else:
            self.blended_img = (self.img_heatmap[:,:,:]*self.blend_slider.get() + self.img_vol*(1 - self.blend_slider.get())).astype('uint8') + 30
            
        img = PIL.Image.fromarray(self.blended_img)
        
        #scaling
        if self.scale_img is not 1.0:
            h,w = img.size
            img = img.resize((int(h*self.scale_img), int(w*self.scale_img)), PIL.Image.ANTIALIAS)
            
        #configure image
        self.photo_blended = PIL.ImageTk.PhotoImage(image = img)
        self.canvasXY.itemconfig(self.image_on_canvas_xy, image = self.photo_blended)
        self.canvasXY.configure(scrollregion = self.canvasXY.bbox("all"))     
        
    def updateVolumeXZ(self):
        """
            Updating XZ Slice by swapping volume
        """   
        self.volumeXZ = np.swapaxes(self.volume, axis1 = 1, axis2 = 0)
        self.updateVolumeXZSlice(None)
    
    def updateVolumeXZSlice(self, event):
        """
            Updating XZ Slice 
        """  
        if self._job_XZ:
            self.master.after_cancel(self._job_XZ)
        self._job_XZ = self.master.after(50, self.updateVolumeXZSlice_Helper)
        
    def updateVolumeXZSlice_Helper(self):
        """
            Updating XZ Slice helper
        """  
        self._job_XZ = None
        #setting up slices
        try:
            self.sliceXZ = self.volumeXZ[self.vol_sliderXZ.get()].astype('uint8')
        except:
            self.sliceXZ = self.volumeXZ[self.volumeXZ.shape[0]//2].astype('uint8')
            
        if(self.SMALLWINDOW == True):
            self.slice_pilXZ = PIL.Image.fromarray(self.sliceXZ[::2,::2])
        else:
            self.slice_pilXZ = PIL.Image.fromarray(self.sliceXZ)
        
        #scaling active
        if self.scale_img is not 1.0:
            h,w = self.slice_pilXZ.size
            self.slice_pilXZ = self.slice_pilXZ.resize((int(h*self.scale_img), int(w*self.scale_img)), PIL.Image.ANTIALIAS)  
        
        #contrast active
        if self.contrast is True:
            self.slice_pilXZ = PIL.ImageEnhance.Contrast(self.slice_pilXZ).enhance(self.contrast_value)
            self.slice_pilXZ = PIL.ImageEnhance.Brightness(self.slice_pilXZ).enhance(self.brightness_value)  
        
        #set up canvas
        try:
            self.photoXZ = PIL.ImageTk.PhotoImage(image = self.slice_pilXZ)
            self.canvasXZ.itemconfig(self.image_on_canvas_XZ, image = self.photoXZ)
            self.canvasXZ.configure(scrollregion = self.canvasXZ.bbox("all")) 
        except:
            return
    
    def updateVolumeYZ(self):
        """
            Updating YZ Slice by swapping volume
        """  
        self.volumeYZ = np.swapaxes(self.volume, axis1 = 0, axis2 = 1)
        self.volumeYZ = np.swapaxes(self.volumeYZ, axis1 = 2, axis2 = 0)
        self.updateVolumeYZSlice(None)
    
    def updateVolumeYZSlice(self, event):
        """
            Updating YZ Slice
            
            Parameters
            ----------
            event: event
                not used but necessary
        """   
        if self._job_YZ:
            self.master.after_cancel(self._job_YZ)
        self._job_YZ = self.master.after(50, self.updateVolumeYZSlice_Helper)
        
    def updateVolumeYZSlice_Helper(self):
        """
            Updating YZ Slice by swapping volume
        """  
        self._job_YZ = None
        #setting up slice
        try:
            self.sliceYZ = self.volumeYZ[self.vol_sliderYZ.get()].astype('uint8')
        except:
            self.sliceYZ = self.volumeYZ[self.volumeYZ.shape[0]//2].astype('uint8')
            
        if(self.SMALLWINDOW == True): 
            self.slice_pilYZ = PIL.Image.fromarray(self.sliceYZ[::2,::2])
        else:
            self.slice_pilYZ = PIL.Image.fromarray(self.sliceYZ)
        #scaling
        if self.scale_img is not 1.0:    
            h,w = self.slice_pilYZ.size
            self.slice_pilYZ = self.slice_pilYZ.resize((int(h*self.scale_img), int(w*self.scale_img)), PIL.Image.ANTIALIAS)    
        #contrast
        if self.contrast is True:
            self.slice_pilYZ = PIL.ImageEnhance.Contrast(self.slice_pilYZ).enhance(self.contrast_value)
            self.slice_pilYZ = PIL.ImageEnhance.Brightness(self.slice_pilYZ).enhance(self.brightness_value) 
        #load image to canvas
        try:    
            self.photoYZ = PIL.ImageTk.PhotoImage(image = self.slice_pilYZ)
            self.canvasYZ.itemconfig(self.image_on_canvas_YZ, image = self.photoYZ)
            self.canvasYZ.configure(scrollregion = self.canvasYZ.bbox("all")) 
        except:
            return
    
    def runRPERefinement(self):
        """ 
            Run RPE refinement
            Extended graph-cut neighborhood (7 neighbors) with limited search space (above BM-only)
            This ensures that no CC is segmented falsely as RPE and improves accuracy.
        """
        if np.count_nonzero(self.segmentation == self.BM_VALUE) == 0:
            print("No BM segmentation present.")
            return
        
        self.statusText.set("Running RPE Refinement...")
        
        if self.shiftedValuesSlowScan is None:
            rfm = RPERefinementAlgorithm.RPERefinementClass(self.volume_original, self.segmentation, self.shiftedValues)
        else:
            rfm = RPERefinementAlgorithm.RPERefinementClass(self.volume_original, self.segmentation, (self.shiftedValues+self.shiftedValuesSlowScan).astype('int32'))
            
        self.segmentation = rfm.runPipeline()
        
        self.heatmap_rpe_bruchs = Heatmap.calculateHeatmap(self.segmentation)
        self.updateSegmentation()
        del rfm
        self.statusText.set("RPE refinement finished")
    
    def runAutomaticSegmentation(self):
        """ 
            Run Automatic three-layer segmentation
            Algorithms can be inserted here - Replace "tls" by own Algorithm and store result in self.segmentation
        """
        self.statusText.set("Running 3-Layer segmentation...")
        
        #INSERT CODE HERE - call segmentation and store
        tls = ThreeLayerSegmentation.LayerSegmentation(self.volume_original, flattening_factor=self.flattening_polynomial,statusText = self.statusText)
        self.segmentation = tls.runPipeline()
        
        # YOU SHOULD NOT TOUCH THE REST HERE 
        #resetting flattening
        self.flatten = True
        self.shiftedValues = np.zeros((self.shiftedValues.shape[0],self.shiftedValues.shape[1]))
        #calculate heatmaps from segmentation
        self.heatmap_rpe_bruchs = Heatmap.calculateHeatmap(self.segmentation)
        self.updateSegmentation()
        del tls
        self.statusText.set("3-Layer Segmentation finished.")
    
    def markGA(self):
        """ 
            Delineate GA selected - Setting up RPE slab
        """
        self.statusText.set("Running 3-Layer segmentation...")
        self.thicknessMarkGA = 2
        self.subBMSlab(self.projMode, int(self.thickness_var.get()))
        self.heatmapType = 'MarkGaByRect' 
        self.heatmapThreshed = self.heatmap.copy()
        self.heatmapForSaving = self.heatmap.copy()
        DelineateGA.pack(self)
    
    def HMDrusen(self): 
        """ 
            Drusen Map (RPEDC) selected
        """
        self.statusText.set("Setting up RPEDC Map...")
        DelineateGA.clearMarkingData(self)
        self.heatmapType = 'Drusen'
        self.updateSegmentation()
        self.statusText.set("RPEDC Map selected")
        
    def minProjectionILMRPE(self):
        """
            Projection ILM to RPE selected
        """   
        self.statusText.set("Setting up ILM to RPE Projection...")
        DelineateGA.clearMarkingData(self)
        self.heatmapType = 'minILMRPEProjection' 
        self.blendSlices()
        self.heatmap = Projections.getProjectionILMRPE(np.asarray(self.volume_original), self.segmentation, self.projMode)
        self.heatmapForSaving = self.heatmap.copy()
        self.heatmap = (255*self.heatmap/65535.0).astype('uint8')
        self.updateSegmentation()
        self.statusText.set("ILM to RPE Projection loaded.")

    def showVessels(self):
        """
            Vesselness selected
        """   
        self.statusText.set("Setting up Vessel image...")
        DelineateGA.clearMarkingData(self)
        self.heatmapType = 'Vesselness'
        self.blendSlices()
        self.vessel_original = Projections.getVesselnessImage(np.asarray(self.volume_original), self.segmentation)
        self.heatmap = self.vessel_original
        self.heatmapForSaving = self.vessel_original.copy()
        self.updateSegmentation()
        self.statusText.set("Vessel image loaded.")
        
    def subRPESlab(self, ga_proj_mode=None, ga_thickness = None): 
        """
            Sub-RPE Slab selected
            
            Optional:
            ---------
            
            ga_proj_mode: string
                projection type(max,min,mean,median)
            ga_thickness: scalar
                thickness of slab
                
        """   
        self.statusText.set("Setting up sub-RPE slab...")
        DelineateGA.clearMarkingData(self)
        self.heatmapType = 'subRPESlab'        
        self.blendSlices()
        # normal and ga delineating mode
        if ga_proj_mode is None and ga_thickness is None:
            try:
                thickness = int(self.thickness_var.get())
            except:
                thickness = 6
                self.thickness_var.set(6)
            self.heatmap = Projections.getsubRPESlab(np.asarray(self.volume_original), self.segmentation, thickness,self.projMode)
        else:
            self.heatmap = Projections.getsubRPESlab(np.asarray(self.volume_original), self.segmentation, ga_thickness, ga_proj_mode)
        #store heatmap
        self.heatmapForSaving = self.heatmap.copy().astype('uint16')
        self.heatmap = (255*self.heatmap/65535.0).astype('uint8')
        self.updateSegmentation()
        self.statusText.set("Sub-RPE slab loaded.")
        
    def subBMSlab(self, ga_proj_mode=None, ga_thickness = None):
        """
            Sub-BM Slab selected
            
            Optional:
            ---------
            
            ga_proj_mode: string
                projection type(max,min,mean,median)
            ga_thickness: scalar
                thickness of slab
                
        """  
        self.statusText.set("Setting up sub-BM slab...")
        self.heatmapType = 'subRPESlab'
        self.blendSlices()
        
        if ga_proj_mode is None and ga_thickness is None:        
            try:
                thickness = int(self.entryThickness.get())
            except:
                thickness = 6
                self.thickness_var.set(6)
            self.heatmap = Projections.getsubBMSlab(np.asarray(self.volume_original), self.segmentation, thickness,self.projMode)
        else:
            self.heatmap = Projections.getsubBMSlab(np.asarray(self.volume_original), self.segmentation, ga_thickness, ga_proj_mode)

        #store heatmap
        self.heatmapForSaving = self.heatmap.copy()
        self.heatmap = (255*self.heatmap/65535.0).astype('uint8')
        self.updateSegmentation()
        self.statusText.set("Sub-BM slab loaded.")
    
    def inpaintOriginalSegmentation(self):
        """
            Inpaint segmentation - Thickness is 3
        """
        self.statusText.set("Inpainting segmentation lines...")

        #flatten original volume (note the minus at shifted volume!!!)
        inpaintedVolume = self.volume_original.copy()
        inpaintedVolume = unFlatten(inpaintedVolume, -self.shiftedValues)

        #volume to rgb
        inpaintedVolume = ImportHandler.getOriginalRGBVolume(inpaintedVolume)
        #segment the three layers
        coordinates_rpe = np.where(self.segmentation == self.RPE_VALUE)
        coordinates_vit = np.where(self.segmentation == self.ILM_VALUE)
        coordinates_bruchs = np.where(self.segmentation == self.BM_VALUE)

        #inpaint with thickness 3
        inpaintedVolume[coordinates_rpe[1],coordinates_rpe[0],coordinates_rpe[2],:] = tuple(int(self.colorRPE.get().lstrip('#')[i:i+2], 16) for i in (0, 2 ,4))
        inpaintedVolume[coordinates_rpe[1]+1,coordinates_rpe[0],coordinates_rpe[2],:] = tuple(int(self.colorRPE.get().lstrip('#')[i:i+2], 16) for i in (0, 2 ,4))
        inpaintedVolume[coordinates_rpe[1]-1,coordinates_rpe[0],coordinates_rpe[2],:] = tuple(int(self.colorRPE.get().lstrip('#')[i:i+2], 16) for i in (0, 2 ,4))
        
        inpaintedVolume[coordinates_bruchs[1],coordinates_bruchs[0],coordinates_bruchs[2],:] = tuple(int(self.colorBM.get().lstrip('#')[i:i+2], 16) for i in (0, 2 ,4))
        inpaintedVolume[coordinates_bruchs[1]+1,coordinates_bruchs[0],coordinates_bruchs[2],:] = tuple(int(self.colorBM.get().lstrip('#')[i:i+2], 16) for i in (0, 2 ,4))
        inpaintedVolume[coordinates_bruchs[1]-1,coordinates_bruchs[0],coordinates_bruchs[2],:] = tuple(int(self.colorBM.get().lstrip('#')[i:i+2], 16) for i in (0, 2 ,4))
        
        inpaintedVolume[coordinates_vit[1],coordinates_vit[0],coordinates_vit[2],:] = tuple(int(self.colorILM.get().lstrip('#')[i:i+2], 16) for i in (0, 2 ,4))
        inpaintedVolume[coordinates_vit[1]+1,coordinates_vit[0],coordinates_vit[2],:] = tuple(int(self.colorILM.get().lstrip('#')[i:i+2], 16) for i in (0, 2 ,4))
        inpaintedVolume[coordinates_vit[1]-1,coordinates_vit[0],coordinates_vit[2],:] = tuple(int(self.colorILM.get().lstrip('#')[i:i+2], 16) for i in (0, 2 ,4))
        
        self.volume = inpaintedVolume
        #update planes
        self.updateVolumeSliceXY(None)
        self.updateVolumeYZ()
        self.updateVolumeXZ()
        self.statusText.set("Inpainting finished.")

    def inpaintSegmentationFiltered(self):
        """
            Inpaint filtered segmentation - Thickness is 3
            When thresholding the heatmap, the results can be backpropagated. 
            This reduces the sensitivity of the map.
        """
        self.statusText.set("Inpainting filtered segmentation...")
        #load original volume
        inpaintedVolume = self.volume_original.copy()
        #flatten original volume (note the minus at shifted volume!!!)
        inpaintedVolume = unFlatten(inpaintedVolume, -self.shiftedValues)
        #volume to rgb
        inpaintedVolume = ImportHandler.getOriginalRGBVolume(inpaintedVolume)
        #inpaint   
        for z in range(self.segmentation.shape[0]):
            coordinates_bruchs = np.where(self.segmentation[z] == self.BM_VALUE)
            
            inpaintedVolume[coordinates_bruchs[0],z,coordinates_bruchs[1],:] = tuple(int(self.colorBM.get().lstrip('#')[i:i+2], 16) for i in (0, 2 ,4))
            inpaintedVolume[coordinates_bruchs[0]+1,z,coordinates_bruchs[1],:] = tuple(int(self.colorBM.get().lstrip('#')[i:i+2], 16) for i in (0, 2 ,4))
            inpaintedVolume[coordinates_bruchs[0]-1,z,coordinates_bruchs[1],:] = tuple(int(self.colorBM.get().lstrip('#')[i:i+2], 16) for i in (0, 2 ,4))
              
            for i in range(coordinates_bruchs[0].size):  
                
                coordinates_rpe_fil = coordinates_bruchs[0][i] - self.heatmapThreshed[z,coordinates_bruchs[1][i]] - 1
            
                inpaintedVolume[coordinates_rpe_fil, z , coordinates_bruchs[1][i],:] =  tuple(int(self.colorRPE.get().lstrip('#')[i:i+2], 16) for i in (0, 2 ,4))
                inpaintedVolume[coordinates_rpe_fil+1, z , coordinates_bruchs[1][i],:] =  tuple(int(self.colorRPE.get().lstrip('#')[i:i+2], 16) for i in (0, 2 ,4))
                inpaintedVolume[coordinates_rpe_fil-1, z , coordinates_bruchs[1][i],:] =  tuple(int(self.colorRPE.get().lstrip('#')[i:i+2], 16) for i in (0, 2 ,4))
        #update planes
        self.volume = inpaintedVolume    
        self.updateVolumeSliceXY(None)
        self.updateVolumeYZ()
        self.updateVolumeXZ()
        self.statusText.set("Inpainting finished.")
        
    def undoSegmentation(self):
        """
            Reset inpainted layers - back to volume only
        """
        self.statusText.set("Reset inpainting...")
        
        inpaintedVolume = self.volume_original.copy()
        #flatten original volume (note the minus at shifted volume!!!)
        inpaintedVolume = unFlatten(inpaintedVolume, -self.shiftedValues)
        #vol to rgb
        inpaintedVolume = ImportHandler.getOriginalRGBVolume(inpaintedVolume)
        #update planes
        self.volume = inpaintedVolume 
        self.updateVolumeSliceXY(None)
        self.updateVolumeXZ()
        self.updateVolumeYZ()
        self.statusText.set("Segmentation reset finished.")
        
    def metricCalculation(self):
        """ 
            Calculate Drusen metrics and show in new window 
        """
        self.metric_map = Heatmap.calculateMetrics(self.heatmap, self.heat_slider_min.get(), self.tv_res.get(), self.axial_res.get())
        self.heatmapType = 'DrusenMetric'
        self.updateSegmentation()
        
    def start_correction_mode(self):
        """ 
            Start correction method 
        """
        #destroy explore windows
        self.FrameVolumeXY.destroy()
        self.FrameVolumeXZ.destroy()
        self.FrameVolumeYZ.destroy()
        self.FrameHeatmap.destroy()
        self.FrameButton.destroy()
        #reset segmentation
        self.segmentation_original = self.segmentation.copy()
        #call manual correction
        try:
            createManualCorrectionFrame(self)
        except Exception as e:
            print(e)