""" 
OptionsPane: Module to create menu panel
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
from PIL import ImageTk, Image, ImageEnhance
import numpy as np
from FileHandler import ImportHandler
from Algorithms.Flattening.VolumeFlattening import applyFlattening,applyFlatteningBM
import threading
from Tooltips.TipHandler import CreateToolTip
import os 

def createOptionMenu(self):
    """ 
    Function to create option menu.
    
    Basically, when new algorithm are inserted, you can call them from here.
    You can either define new Buttons/RadioButtons or substitute the methods
    for the old Buttons.
    """
    def background(func):
        """ 
            background worker funcrion - when loading function, gui doesnt stop
            
            Parameters
            ----------
            func: function
                function
        """
        th = threading.Thread(target=func)
        th.start()
    
    def createStatusLabel():
        """ 
           Status label to track progess
        """        
        self.statusText = StringVar()
        self.statusText.set("")
        self.LabelPGBar = Label(self.FrameButton, width = self.btn_double_width,height= 2, textvariable=self.statusText, font= self.statusFont,anchor="n", bg=self.canvasbackground, fg =self.hlghtbg)
        self.LabelPGBar.grid(row=0,column=0,columnspan=2,pady=(25,0))
    
    def applyFlatteningBM_Helper():
        """ 
            Flattening to BM
        """ 
        self.statusText.set("Flattening to BM...")
        applyFlatteningBM(self)
        self.statusText.set("BM flattening done!")
        
    def applyFlattening_Helper():
        """ 
            Flattening to RPE
        """ 
        self.statusText.set("Flattening to RPE...")
        applyFlattening(self)
        self.statusText.set("RPE flattening done!")

    def threshToOne():
        """ 
            Button callback to thresh to black/white or black/gray
        """ 
        if self.threshold_to_one is True:
            self.ButtonThreshToOne.configure(text="Threshold To White")
            self.threshold_to_one = False
        else:
            self.ButtonThreshToOne.configure(text="Threshold To Gray")
            self.threshold_to_one = True
    
    def setThreshold():
        """ 
            Apply thresholding with value set in self.scale_threshold
        """ 
        if self.threshold_activate is True:
            self.undoSegmentation()
            self.statusText.set("Thresholding...")
            if self.threshold_to_one is False:
                self.volume = ImportHandler.getOriginalRGBVolume(np.where(self.volume_original < int(self.scale_threshold.get()),0,self.volume_original))
            else:
                maxval = np.max(self.volume_original) 
                self.volume = ImportHandler.getOriginalRGBVolume(np.where(self.volume_original < int(self.scale_threshold.get()),0, maxval))
            self.statusText.set("Updating planes...")
            self.updateVolumeSliceXY(None)
            self.updateVolumeXZ()
            self.updateVolumeYZ()
            self.threshold_activate = False
            self.ButtonThresholdActivate.configure(text="Deactivate", bg=self.btn_common_bg)
            self.statusText.set("Thresholding done!")
        else:
            self.statusText.set("Reset Thresholding...")
            self.volume = ImportHandler.getOriginalRGBVolume(self.volume_original)
            self.statusText.set("Updating planes...")
            self.updateVolumeSliceXY(None)
            self.updateVolumeXZ()
            self.updateVolumeYZ()
            self.threshold_activate = True
            self.ButtonThresholdActivate.configure(text="Activate", bg=self.btnbackground)
            self.statusText.set("Threshold reset done!")
    
    def colormapChange(selection):
        """ 
            Change colormap of heatmap
            
            Parameters
            ----------
            selection: event, StringVar
                selected layer from dropdown menu
        """ 
        if selection is 'Viridis':
            self.COLORMAP = 'viridis'
        elif selection is 'Diverging':
            self.COLORMAP = 'RdBu_r'
        elif selection is 'Jet':
            self.COLORMAP = 'jet'
        elif selection is 'Grayvalues':
            self.COLORMAP = 'gray'
        self.statusText.set("Colormap: "+selection)
        self.updateHeatSlice(None)
        
    def colorLineChanged(event=None):
        """ 
            Segmentation Line Color changed and updated.
            
            Optional
            ----------
            event: event
                unused event
        """ 
        self.statusText.set("Updating Line colors")
        background(self.inpaintOriginalSegmentation)
        self.canvasColorBM.configure(bg=self.colorBM.get())
        self.canvasColorRPE.configure(bg=self.colorRPE.get())
        self.canvasColorILM.configure(bg=self.colorILM.get())
        
    def projectionModeChange(selection):
        """ 
            Projection mode selection
            
            Parameters
            ----------
            selection: event, StringVar
                selected layer from dropdown menu
        """ 
        if selection is 'Mean':
            self.projMode = 'mean'
        elif selection is 'Median':
            self.projMode = 'median'
        elif selection is 'Maximum':
            self.projMode = 'max'
        elif selection is 'Minimum':
            self.projMode = 'min'
        self.statusText.set("Projection mode: "+selection)
        
        if(self.varVis.get() == 3):
            background(self.subRPESlab)
        elif(self.varVis.get() == 4):
            background(self.subBMSlab)
        elif(self.varVis.get() == 6):
            background(self.minProjectionILMRPE)
     
    def thicknessChanged(event=None):
        """
            Callback when thickness changed
            
            Optional
            ----------
            event: event
                unused event
        """
        self.statusText.set("Thickness updated")
        
        if(self.varVis.get() == 3):
            background(self.subRPESlab)
        elif(self.varVis.get() == 4):
            background(self.subBMSlab)
        elif(self.varVis.get() == 6):
            background(self.minProjectionILMRPE)   
    
    def slideContrast(event = None):
        """ 
            Slider callback to set contrast
            
            Optional
            ----------
            event: event
                unused event
        """ 
        if self.contrast is True:
            self.contrast_value = self.contrast_slider.get()
            self.brightness_value = self.brightness_slider.get()
            self.updateVolumeSliceXY(None)
            self.updateVolumeYZSlice(None)
            self.updateVolumeXZSlice(None)    
    
    def activateContrastMenu():
        """ 
            (De-)Activate contrast setting
        """ 
        if self.contrast is True:
            self.statusText.set("Resetting contrast...")
            self.buttonContrast.configure(text="Activate", bg=self.btnbackground)
            self.contrast_value = self.contrast_slider.get()
            self.brightness_value = self.brightness_slider.get()
            self.updateVolumeSliceXY(None)
            self.updateVolumeYZSlice(None)
            self.updateVolumeXZSlice(None)    
            self.contrast = False
            self.statusText.set("Contrast reset done!")
        else:
            self.statusText.set("Contrast variation active!")
            self.buttonContrast.configure(text="Deactivate", bg=self.btn_common_bg)
            self.contrast_value = self.contrast_slider.get()
            self.brightness_value = self.brightness_slider.get()
            self.updateVolumeSliceXY(None)
            self.updateVolumeYZSlice(None)
            self.updateVolumeXZSlice(None)  
            self.contrast = True
        
    def blendButton():
        """ 
            (De-)Activate blending to en face plane
        """ 
        if self.blender is True:
            self.statusText.set("Blender turned off.")
            self.blenderButton.configure(text="Activate", bg=self.btnbackground)
            self.updateVolumeSliceXY(None)
            self.blender = False
        else:
            self.statusText.set("Blender active!")
            self.blenderButton.configure(text="Deactivate", bg=self.btn_common_bg)
            self.blender = True
            self.blendSlices() 
        
    def updateMetrics():
        """ 
            Update metrics when changing resolution
        """ 
        self.statusText.set("Calculating metrics...")
        try:
            val = np.sum(np.where(self.heatmapThreshed != 0, 1, 0))
            self.labelDrusenArea.configure(text=str(np.round(val*self.tv_res.get()**2*1e-6,8)))
            self.labelDrusenAreaPercent.configure(text=str(np.round(val/(self.heatmapThreshed.shape[0]*self.heatmapThreshed.shape[1]),8)*100))
            
            val = np.sum(self.heatmapThreshed)
            self.labelDrusenVolume.configure(text= str(np.round(val*self.tv_res.get()**2*self.axial_res.get()*1e-9,8)))
            self.labelDrusenVolumePercent.configure(text= str(np.round(val/(self.volume.shape[0]*self.volume.shape[1]*self.volume.shape[2]),8)*100))
        except:
            pass
        self.statusText.set("")    
    
    def callDataMenu():
        """ 
            Callback when data button clicked.
            
            Constructing Data Menu.
        """ 
        if self.FrameButton.winfo_exists():
            self.FrameButton.destroy()
            
        self.FrameButton = Frame(self.master, width = self.menuwidth , height = self.option_height, bg=self.canvasbackground)
        self.FrameButton.grid(row = 0, column = 1, rowspan=2, sticky='nw')
    
        createStatusLabel()
        
        Label(self.FrameButton, width = self.btn_single_width,height= 1, text="OCT Volume",font=self.headerFont,anchor="n", bg=self.canvasbackground, fg =self.canvasforeground).grid(row=0,column=0,columnspan=2,pady=self.frame_pad_y)
        self.btnLoadUnmerged = Button(self.FrameButton,width = self.btn_single_width , height = self.btn_height , text="Load w/ Flattening", font=self.btnFont, command=lambda:background(self.loadVolumeUnmerged), fg=self.btnforeground, bg=self.btnbackground, activebackground=self.btn_common_bg, activeforeground=self.btnforeground)
        self.btnLoadUnmerged.grid(row=1,column=0,padx=self.pad_x_left,pady=(10,0), sticky="W")
        CreateToolTip(self.btnLoadUnmerged, self.ttip_dict['loadRaw'])
        
        self.btnLoadMerged = Button(self.FrameButton,width = self.btn_single_width , height = self.btn_height , text="Load w/o Flattening", font=self.btnFont, command=lambda:background(self.loadVolume), fg=self.btnforeground, bg=self.btnbackground, activebackground=self.btn_common_bg, activeforeground=self.btnforeground)
        self.btnLoadMerged.grid(row=2,column=0,padx=self.pad_x_left,pady=(10,0), sticky="W")
        CreateToolTip(self.btnLoadMerged, self.ttip_dict['loadMerged'])
        
        self.btnLoadOCTA = Button(self.FrameButton,width = self.btn_single_width , height = self.btn_height , text="Load OCTA", font=self.btnFont, command=lambda:background(self.loadVolumeOCTA), fg=self.btnforeground, bg=self.btnbackground, activebackground=self.btn_common_bg, activeforeground=self.btnforeground)
        self.btnLoadOCTA.grid(row=1,column=1,padx=(5,0),pady=(10,0), sticky="W")
        CreateToolTip(self.btnLoadOCTA, self.ttip_dict['loadOCTA'])
        
        Label(self.FrameButton, width = self.btn_single_width,height= 1, text="Segmentation",font=self.headerFont,anchor="n", bg=self.canvasbackground, fg =self.canvasforeground).grid(row=3,column=0,columnspan=2, pady=self.separate_label_pad_y)
        self.btnLoadSeg = Button(self.FrameButton,width = self.btn_single_width , height = self.btn_height , text="Load Data", font=self.btnFont, command=lambda:background(self.loadCompleteSegmentation), fg=self.btnforeground, bg=self.btnbackground, activebackground=self.btn_common_bg, activeforeground=self.btnforeground)
        self.btnLoadSeg.grid(row=4,column=0,padx=self.pad_x_left,pady=(10,0), sticky="W")
        CreateToolTip(self.btnLoadSeg, self.ttip_dict['loadSeg'])
        
        self.btnResetSeg = Button(self.FrameButton,width = self.btn_single_width , height = self.btn_height , text="Reset Data", font=self.btnFont, command=lambda:background(self.resetSegmentation), fg=self.btnforeground, bg=self.btnbackground, activebackground=self.btn_common_bg, activeforeground=self.btnforeground)
        self.btnResetSeg.grid(row=4,column=1,padx=(5,0),pady=(10,0), sticky="W")
        CreateToolTip(self.btnResetSeg, self.ttip_dict['resetSeg'])
        
        self.btnSaveSegInpainted = Button(self.FrameButton,width = self.btn_single_width , height = self.btn_height , text="Save Inpainted", font=self.btnFont, command=lambda:background(self.saveInpainted), fg=self.btnforeground, bg=self.btnbackground, activebackground=self.btn_common_bg, activeforeground=self.btnforeground)
        self.btnSaveSegInpainted.grid(row=5,column=0,padx=self.pad_x_left,pady=(10,0), sticky="W")
        CreateToolTip(self.btnSaveSegInpainted, self.ttip_dict['saveInpainted'])
        
        self.btnSaveSegLines = Button(self.FrameButton,width = self.btn_single_width , height = self.btn_height , text="Save Lines", font=self.btnFont, command=lambda:background(self.saveLineSegmentation), fg=self.btnforeground, bg=self.btnbackground, activebackground=self.btn_common_bg, activeforeground=self.btnforeground)
        self.btnSaveSegLines.grid(row=5,column=1,padx=(5,0),pady=(10,0), sticky="W")
        CreateToolTip(self.btnSaveSegLines, self.ttip_dict['saveLines'])
        
        Label(self.FrameButton, width = self.btn_single_width,height= 1, text="Visualization",font=self.headerFont,anchor="n", bg=self.canvasbackground, fg =self.canvasforeground).grid(row=6,column=0,columnspan=2, pady=self.separate_label_pad_y)
        self.btnSaveHM = Button(self.FrameButton,width = self.btn_double_width , height = self.btn_height , text="Save Image", font=self.btnFont, command=lambda:background(self.saveHeatmap), fg=self.btnforeground, bg=self.btnbackground, activebackground=self.btn_common_bg, activeforeground=self.btnforeground)
        self.btnSaveHM .grid(row=7,column=0, columnspan=2,padx=self.pad_x_right,pady=(10,0), sticky="W")
        CreateToolTip(self.btnSaveHM, self.ttip_dict['saveHM'])

    def callImageMenu():
        """ 
            Callback when Image enhancement button clicked.
            
            Constructing Image Enhancement Menu.
        """ 
        if self.FrameButton.winfo_exists():
            self.FrameButton.destroy()
        
        self.FrameButton = Frame(self.master, width = self.menuwidth , height = self.option_height, bg=self.canvasbackground)
        self.FrameButton.grid(row = 0, column = 1, rowspan=2, sticky='nw')
        
        createStatusLabel()
        
        Label(self.FrameButton, width = self.btn_single_width,height= 1, text="Image Enhancement",font=self.headerFont,anchor="n", bg=self.canvasbackground, fg =self.canvasforeground).grid(row=0,column=0,columnspan=2,pady=self.frame_pad_y)
        
        Label(self.FrameButton, width = 20,height= 1, text="Brightness",font=self.textFont,anchor="n", bg=self.canvasbackground, fg =self.canvasforeground).grid(row=1,column=0,columnspan=2,pady=self.separate_label_pad_y)
        self.brightness_slider = Scale(self.FrameButton, from_=0, to=5, resolution=0.05, orient=HORIZONTAL, command=slideContrast,length = self.scale_length, font = self.textFont, bg=self.canvasbackground, fg =self.canvasforeground )
        self.brightness_slider.grid(row=2,column=0, columnspan = 2, sticky =NW, padx=self.pad_x_left,pady=(10,0))
        self.brightness_slider.set(1.0)
        
        Label(self.FrameButton, width = 20,height= 1, text="Contrast",font=self.textFont,anchor="n", bg=self.canvasbackground, fg =self.canvasforeground).grid(row=3,column=0,columnspan=2,pady=self.separate_label_pad_y)
        self.contrast_slider = Scale(self.FrameButton, from_=0, to=5, resolution=0.05, orient=HORIZONTAL, command=slideContrast,length = self.scale_length, font = self.textFont, bg=self.canvasbackground, fg =self.canvasforeground )
        self.contrast_slider.grid(row=4,column=0, columnspan = 2, sticky =NW, padx=self.pad_x_left,pady=(10,0))
        self.contrast_slider.set(1.0)
        
        if(self.contrast is False):
            self.buttonContrast = Button(self.FrameButton,width = self.btn_double_width, height = self.btn_height , text="Activate", font=self.btnFont, command=activateContrastMenu, fg=self.btnforeground, bg=self.btnbackground)
        else:
            self.buttonContrast = Button(self.FrameButton,width = self.btn_double_width, height = self.btn_height , text="Deactivate", font=self.btnFont, command=activateContrastMenu, fg=self.btnforeground,bg=self.btn_common_bg)
        self.buttonContrast.grid(row=5,column=0,columnspan=2,padx=self.pad_x_right,pady=(10,0), sticky="W")
        CreateToolTip(self.buttonContrast, self.ttip_dict['contrast'])
        
        Label(self.FrameButton, width = self.btn_single_width,height= 1, text="Thresholding",font=self.headerFont,anchor="n", bg=self.canvasbackground, fg =self.canvasforeground).grid(row=6,column=0,columnspan=2,pady=self.separate_label_pad_y)
        self.scale_threshold = Scale(self.FrameButton, from_=0, to=np.max(self.volume_original), length = self.scale_length,resolution=1.0, orient=HORIZONTAL, font = self.textFont, bg=self.canvasbackground, fg =self.canvasforeground)
        self.scale_threshold.grid(row=7,column=0,columnspan=2,padx=self.pad_x_left,pady=(10,0),sticky="W")
        
        self.ButtonThreshToOne = Button(self.FrameButton,width = self.btn_double_width , height = self.btn_height , text="Threshold To White", font=self.btnFont, command=threshToOne, fg=self.btnforeground, bg=self.btnbackground)
        self.ButtonThreshToOne.grid(row=8,column=0,columnspan=2,padx=self.pad_x_right,pady=(10,0),sticky="W")
        CreateToolTip(self.ButtonThreshToOne,self.ttip_dict['thresholdToOne'])
        if(self.threshold_activate is True):
            self.ButtonThresholdActivate = Button(self.FrameButton,width = self.btn_double_width , height = self.btn_height , text="Activate", font=self.btnFont, command=lambda: background(setThreshold), fg=self.btnforeground, bg=self.btnbackground)
        else:
            self.ButtonThresholdActivate = Button(self.FrameButton,width = self.btn_double_width , height = self.btn_height , text="Deactivate", font=self.btnFont, command=setThreshold, fg=self.btnforeground, bg=self.btn_common_bg)
        self.ButtonThresholdActivate.grid(row=9,column=0,columnspan=2,padx=self.pad_x_right,pady=(10,0), sticky="W")
        CreateToolTip(self.ButtonThresholdActivate, self.ttip_dict['threshold'])
    
        self.varColormap = StringVar(self.FrameButton)
        self.varColormap.set("Viridis")
        Label(self.FrameButton, width = self.btn_single_width,height= 1, text="Colormap",font=self.headerFont,anchor="n", bg=self.canvasbackground, fg =self.canvasforeground).grid(row=10,column=0,columnspan=2,pady=self.separate_label_pad_y)
        self.optionColormap = OptionMenu(self.FrameButton, self.varColormap, "Viridis",  "Diverging", "Jet", "Grayvalues", command=colormapChange)
        self.optionColormap.config(font=self.textFont)
        menu = self.optionColormap.nametowidget(self.optionColormap.menuname)
        menu.configure(font=self.textFont)
        self.optionColormap.grid(row=11,column=0, columnspan=2,pady=(10,0), sticky="N")
        CreateToolTip(self.optionColormap, self.ttip_dict['colormap'])
        
    def callSegmentationMenu():
        """ 
            Callback when Segmentation button clicked.
            
            Constructing Segmentation Menu.
        """         
        if self.FrameButton.winfo_exists():
            self.FrameButton.destroy()
        
        ctr_rows = 0
        
        self.FrameButton = Frame(self.master, width = self.menuwidth , height = self.option_height, bg=self.canvasbackground)
        self.FrameButton.grid(row = 0, column = 1, rowspan=2, sticky='nw')
                    
        createStatusLabel()
        
        Label(self.FrameButton, width = self.btn_single_width,height= 1, text="GUI Mode",font=self.headerFont,anchor="n", bg=self.canvasbackground, fg =self.canvasforeground).grid(row=ctr_rows, column=0, columnspan=2, pady=self.frame_pad_y)
        ctr_rows += 1
        self.btnManRefine = Button(self.FrameButton,width = self.btn_double_width , height = self.btn_height , text="Manual Refinement", font=self.btnFont, command=self.correctionmode, fg=self.btnforeground, bg=self.btnbackground, activebackground=self.btn_common_bg, activeforeground=self.btnforeground)
        self.btnManRefine.grid(row=ctr_rows,column=0,columnspan=2,padx=self.pad_x_right,pady=self.separate_label_pad_y, sticky="W")
        CreateToolTip(self.btnManRefine, self.ttip_dict['manRefine'])
        ctr_rows += 1
        Label(self.FrameButton, width = 20,height= 1, text="Automatic Segmentation",font=self.headerFont,anchor="n", bg=self.canvasbackground, fg =self.canvasforeground).grid(row=ctr_rows, column=0, columnspan=2, pady=self.separate_label_pad_y)
        ctr_rows += 1
        self.btnThreeLayer = Button(self.FrameButton,width = self.btn_single_width , height = self.btn_height , text="3-layer Algorithm", font=self.btnFont, command=lambda:background(self.runAutomaticSegmentation), fg=self.btnforeground, bg=self.btnbackground, activebackground=self.btn_common_bg, activeforeground=self.btnforeground)
        self.btnThreeLayer.grid(row=ctr_rows,column=0,padx=self.pad_x_left,pady=self.separate_label_pad_y, sticky="W")
        CreateToolTip(self.btnThreeLayer, self.ttip_dict['threeLayer'])
        self.btnRPERefine = Button(self.FrameButton,width = self.btn_single_width , height = self.btn_height , text="RPE Refinement", font=self.btnFont, command=lambda:background(self.runRPERefinement), fg=self.btnforeground, bg=self.btnbackground, activebackground=self.btn_common_bg, activeforeground=self.btnforeground)
        self.btnRPERefine.grid(row=ctr_rows,column=1,padx=(5,0),pady=self.separate_label_pad_y, sticky="W")
        CreateToolTip(self.btnRPERefine, self.ttip_dict['RPERefine'])
        ctr_rows += 1
        
        Label(self.FrameButton, width = 20,height= 1, text="Inpainting Segmentation",font=self.headerFont,anchor="n", bg=self.canvasbackground, fg =self.canvasforeground).grid(row=ctr_rows, column=0, columnspan=2, pady=self.separate_label_pad_y)   
        ctr_rows += 1     
        self.btnInpaintOrg = Button(self.FrameButton,width = self.btn_single_width , height = self.btn_height , text="Inpaint Original", font=self.btnFont, command=lambda:background(self.inpaintOriginalSegmentation), fg=self.btnforeground, bg=self.btnbackground, activebackground=self.btn_common_bg, activeforeground=self.btnforeground)
        self.btnInpaintOrg.grid(row=ctr_rows,column=0,padx=self.pad_x_left,pady=self.separate_label_pad_y, sticky="W")
        CreateToolTip(self.btnInpaintOrg, self.ttip_dict['segInpaint'])
        self.InpaintFil = Button(self.FrameButton,width = self.btn_single_width , height = self.btn_height , text="Inpaint Filtered", font=self.btnFont, command=lambda:background(self.inpaintSegmentationFiltered), fg=self.btnforeground, bg=self.btnbackground, activebackground=self.btn_common_bg, activeforeground=self.btnforeground)
        self.InpaintFil.grid(row=ctr_rows,column=1,padx=(5,0),pady=self.separate_label_pad_y, sticky="W")
        CreateToolTip(self.InpaintFil, self.ttip_dict['filInpaint'])
        ctr_rows += 1
        self.InpaintReset = Button(self.FrameButton,width = self.btn_double_width , height = self.btn_height , text="Reset Inpainting", font=self.btnFont, command=lambda:background(self.undoSegmentation), fg=self.btnforeground, bg=self.btnbackground, activebackground=self.btn_common_bg, activeforeground=self.btnforeground)
        self.InpaintReset.grid(row=ctr_rows,column=0, columnspan=2,padx=self.pad_x_right,pady=(10,0), sticky="W")
        CreateToolTip(self.InpaintReset, self.ttip_dict['resetInpaint'])
        ctr_rows += 1
    
        Label(self.FrameButton, width = 20,height= 1, text="Color ILM:",font=self.textFont,anchor="nw", bg=self.canvasbackground, fg =self.canvasforeground).grid(row=ctr_rows,column=0, columnspan=2, padx=self.pad_x_left,pady=(5,0), sticky="W")  
        self.entryColorILM = Entry(self.FrameButton, textvariable=self.colorILM,font=self.textFont)
        self.entryColorILM.grid(row=ctr_rows,column=0, columnspan=2, padx = self.entry_pad_x,pady=(5,0), sticky="NW")
        self.canvasColorILM = Canvas(self.FrameButton, width = 5, height = 5, bg=self.colorILM.get(),highlightthickness = 0)
        self.canvasColorILM.grid(row=ctr_rows,column=0, columnspan=2, padx = (self.entry_pad_x[0]-5,0),pady=(5,0), sticky="NW")
        CreateToolTip(self.entryColorILM, self.ttip_dict['colorILM'])
        self.entryColorILM.bind('<Return>', colorLineChanged)
        ctr_rows += 1
        
        Label(self.FrameButton, width = 20,height= 1, text="Color RPE:",font=self.textFont,anchor="nw", bg=self.canvasbackground, fg =self.canvasforeground).grid(row=ctr_rows,column=0, columnspan=2, padx=self.pad_x_left,pady=(5,0), sticky="W")  
        self.entryColorRPE = Entry(self.FrameButton, textvariable=self.colorRPE,font=self.textFont)
        self.entryColorRPE.grid(row=ctr_rows,column=0, columnspan=2, padx = self.entry_pad_x,pady=(5,0), sticky="NW")
        self.canvasColorRPE = Canvas(self.FrameButton, width = 5, height = 5, bg=self.colorRPE.get(),highlightthickness = 0)
        self.canvasColorRPE.grid(row=ctr_rows,column=0, columnspan=2, padx = (self.entry_pad_x[0]-5,0),pady=(5,0), sticky="NW")
        CreateToolTip(self.entryColorRPE, self.ttip_dict['colorRPE'])
        self.entryColorRPE.bind('<Return>', colorLineChanged)
        ctr_rows += 1

        Label(self.FrameButton, width = 20,height= 1, text="Color BM:",font=self.textFont,anchor="nw", bg=self.canvasbackground, fg =self.canvasforeground).grid(row=ctr_rows,column=0, columnspan=2, padx=self.pad_x_left,pady=(5,0), sticky="W")  
        self.entryColorBM = Entry(self.FrameButton, textvariable=self.colorBM,font=self.textFont)
        self.entryColorBM.grid(row=ctr_rows,column=0, columnspan=2, padx = self.entry_pad_x,pady=(5,0), sticky="NW")
        self.canvasColorBM = Canvas(self.FrameButton, width = 5, height = 5, bg=self.colorBM.get(),highlightthickness = 0)
        self.canvasColorBM.grid(row=ctr_rows,column=0, columnspan=2, padx = (self.entry_pad_x[0]-5,0),pady=(5,0), sticky="NW")
        CreateToolTip(self.entryColorBM, self.ttip_dict['colorBM'])
        self.entryColorBM.bind('<Return>', colorLineChanged)
        ctr_rows += 1
        
        Label(self.FrameButton, width = self.btn_single_width,height= 1, text="Flattening",font=self.headerFont,anchor="n", bg=self.canvasbackground, fg =self.canvasforeground).grid(row=ctr_rows, column=0, columnspan=2, pady=self.separate_label_pad_y) 
        ctr_rows += 1
        
        if(self.flattenBM is True):
            self.buttonFlatteningBM = Button(self.FrameButton,width = self.btn_double_width , height = self.btn_height , text="Flatten to Bruch's Membrane", font=self.btnFont, command=lambda: background(applyFlatteningBM_Helper), fg=self.btnforeground, bg=self.btnbackground, activebackground=self.btn_common_bg, activeforeground=self.btnforeground)
        else:
            self.buttonFlatteningBM = Button(self.FrameButton,width = self.btn_double_width , height = self.btn_height , text="Reset Bruch's Flattening", font=self.btnFont, command=lambda: background(applyFlatteningBM_Helper), fg=self.btnforeground, bg=self.btn_common_bg, activebackground=self.btn_common_bg, activeforeground=self.btnforeground)
        self.buttonFlatteningBM.grid(row=ctr_rows,column=0, columnspan=2,padx=self.pad_x_right,pady=(10,0), sticky="W")
        CreateToolTip(self.buttonFlatteningBM, self.ttip_dict['flattBM'])
        ctr_rows += 1
        
        self.buttonFlattening = Button(self.FrameButton,width = self.btn_double_width , height = self.btn_height , text="Apply Flattening to RPE", font=self.btnFont, command=lambda: background(applyFlattening_Helper), fg=self.btnforeground, bg=self.btnbackground, activebackground=self.btn_common_bg, activeforeground=self.btnforeground)
        self.buttonFlattening.grid(row=ctr_rows,column=0, columnspan=2,padx=self.pad_x_right,pady=(10,0), sticky="W")
        CreateToolTip(self.buttonFlattening, self.ttip_dict['flattRPE'])
        ctr_rows += 1
        
        self.polynomial_var = IntVar()
        self.polynomial_var.set(4)
        Label(self.FrameButton, width = 20,height= 1, text="Order of Polynomial:",font=self.textFont,anchor="nw", bg=self.canvasbackground, fg =self.canvasforeground).grid(row=ctr_rows,column=0, columnspan=2, padx=self.pad_x_left,pady=(10,0), sticky="W")  
        self.entryPolynomial = Entry(self.FrameButton, textvariable=self.polynomial_var, font=self.textFont)
        self.entryPolynomial.grid(row=ctr_rows,column=0, columnspan=2,padx = self.entry_pad_x,pady=(10,0), sticky="NW")

    def callVisualizationMenu():
        """ 
            Callback when Visualization button clicked.
            
            Constructing Visualization Menu.
        """           
        if self.FrameButton.winfo_exists():
            self.FrameButton.destroy()
        
        self.FrameButton = Frame(self.master, width = self.menuwidth , height = self.option_height, bg=self.canvasbackground)
        self.FrameButton.grid(row = 0, column = 1, rowspan=2, sticky='nw')
        
        createStatusLabel()
        
        Label(self.FrameButton, width = self.btn_single_width,height= 1, text="Drusen",font=self.headerFont,anchor="n", bg=self.canvasbackground, fg =self.canvasforeground).grid(row=0, column=0, columnspan=2, pady=self.frame_pad_y)

        self.btnRPEDC = Radiobutton(self.FrameButton,width = self.btn_double_width , height = self.btn_height , text="RPEDC Thickness Map", font=self.btnFont,variable=self.varVis, value=1,indicatoron=0, command=lambda:background(self.HMDrusen), activebackground=self.btn_common_bg, fg=self.btnforeground, bg=self.btnbackground, activeforeground=self.btnforeground, selectcolor=self.btn_common_bg)
        self.btnRPEDC.grid(row=1,column=0, columnspan=2,padx=self.pad_x_right,pady=(10,0), sticky="N")
        CreateToolTip(self.btnRPEDC, self.ttip_dict['rpedcthick'])
        self.btnSevere = Radiobutton(self.FrameButton,width = self.btn_double_width , height = self.btn_height , text="Drusen Severeness", font=self.btnFont,variable=self.varVis, value=2,indicatoron=0, command=self.metricCalculation, fg=self.btnforeground, bg=self.btnbackground, activebackground=self.btn_common_bg, activeforeground=self.btnforeground, selectcolor=self.btn_common_bg)
        self.btnSevere.grid(row=2,column=0, columnspan=2,padx=self.pad_x_right,pady=(10,0), sticky="N")
        CreateToolTip(self.btnSevere, self.ttip_dict['severeness'])

        Label(self.FrameButton, width = self.btn_single_width,height= 1, text="Projections",font=self.headerFont,anchor="n", bg=self.canvasbackground, fg =self.canvasforeground).grid(row=3, column=0, columnspan=2, pady=self.separate_label_pad_y)        
        self.btnsubRPE = Radiobutton(self.FrameButton,width = self.btn_double_width , height = self.btn_height , text="Sub-RPE Slab", font=self.btnFont,variable=self.varVis, value=3,indicatoron=0, command=lambda:background(self.subRPESlab), fg=self.btnforeground, bg=self.btnbackground, activebackground=self.btn_common_bg, activeforeground=self.btnforeground, selectcolor=self.btn_common_bg)
        self.btnsubRPE.grid(row=4,column=0, columnspan=2,padx=self.pad_x_right,pady=(10,0), sticky="N")
        CreateToolTip(self.btnsubRPE, self.ttip_dict['subrpe'])
        self.btnSubBM = Radiobutton(self.FrameButton,width = self.btn_double_width , height = self.btn_height , text="Sub-BM Slab", font=self.btnFont,variable=self.varVis, value=4,indicatoron=0, command=lambda:background(self.subBMSlab), fg=self.btnforeground, bg=self.btnbackground, activebackground=self.btn_common_bg, activeforeground=self.btnforeground, selectcolor=self.btn_common_bg)
        self.btnSubBM.grid(row=5,column=0, columnspan=2,padx=self.pad_x_right,pady=(10,0), sticky="N") 
        CreateToolTip(self.btnSubBM, self.ttip_dict['subbm'])
        self.btnILMRPE = Radiobutton(self.FrameButton,width = self.btn_double_width , height = self.btn_height , text="ILM - RPE", font=self.btnFont,variable=self.varVis, value=6,indicatoron=0, command=lambda:background(self.minProjectionILMRPE), fg=self.btnforeground, bg=self.btnbackground, activebackground=self.btn_common_bg, activeforeground=self.btnforeground, selectcolor=self.btn_common_bg)
        self.btnILMRPE.grid(row=6,column=0, columnspan=2,padx=self.pad_x_right,pady=(10,0), sticky="N")
        CreateToolTip(self.btnILMRPE, self.ttip_dict['ilmrpe'])
        
        self.thickness_var = IntVar()
        self.thickness_var.set(6)
        Label(self.FrameButton, width = 20,height= 1, text="Thickness [px]:",font=self.textFont,anchor="nw", bg=self.canvasbackground, fg =self.canvasforeground).grid(row=7,column=0, columnspan=2, padx=self.pad_x_left,pady=(5,0), sticky="W")  
        self.entryThickness = Entry(self.FrameButton, textvariable=self.thickness_var,font=self.textFont)
        self.entryThickness.grid(row=7,column=0, columnspan=2,padx = self.entry_pad_x,pady=(5,0), sticky="NW")
        CreateToolTip(self.entryThickness, self.ttip_dict['slabthickness'])
        self.entryThickness.bind('<Return>', thicknessChanged)
        
        Label(self.FrameButton, width = 5,height= 1, text="Mode:",font=self.textFont,anchor="n", bg=self.canvasbackground, fg =self.canvasforeground).grid(row=8,column=0,columnspan=2,  padx=self.pad_x_left,pady=(5,0), sticky="W")
        self.varOptionProj = StringVar(self.FrameButton)
        self.varOptionProj.set("Mean")
        self.optionProj = OptionMenu(self.FrameButton, self.varOptionProj, "Mean", "Median", "Minimum", "Maximum", command=projectionModeChange)
        self.optionProj.config(font=self.textFont)
        menu = self.optionProj.nametowidget(self.optionProj.menuname)
        menu.configure(font=self.textFont)
        self.optionProj.grid(row=8,column=0, columnspan=2, padx=(207,0),pady=(5,0), sticky="W")
        CreateToolTip(self.optionProj, self.ttip_dict['projectionmode'])
        
        Label(self.FrameButton, width = self.btn_single_width,height= 1, text="Miscellaneous",font=self.headerFont,anchor="n", bg=self.canvasbackground, fg =self.canvasforeground).grid(row=9, column=0, columnspan=2, pady=self.separate_label_pad_y)       
        self.btnVessel = Radiobutton(self.FrameButton,width = self.btn_double_width , height = self.btn_height , text="Show Vessels", font=self.btnFont,variable=self.varVis, value=7,indicatoron=0, command=lambda:background(self.showVessels), activebackground=self.btn_common_bg, fg=self.btnforeground, bg=self.btnbackground, activeforeground=self.btnforeground, selectcolor=self.btn_common_bg)
        self.btnVessel.grid(row=10,column=0, columnspan=2,padx=self.pad_x_right,pady=(10,0), sticky="N")
        CreateToolTip(self.btnVessel, self.ttip_dict['vesselness'])
        self.btnGA = Radiobutton(self.FrameButton,width = self.btn_double_width , height = self.btn_height , text="Delineate GA", font=self.btnFont, variable=self.varVis, value=8,indicatoron=0, command=self.markGA, fg=self.btnforeground, bg=self.btnbackground, activebackground=self.btn_common_bg, activeforeground=self.btnforeground,selectcolor=self.btn_common_bg)
        self.btnGA.grid(row=11,column=0, columnspan=2,padx=self.pad_x_right,pady=(10,0), sticky="N")
        CreateToolTip(self.btnGA, self.ttip_dict['delineateGA'])
        
        Label(self.FrameButton, width = self.btn_single_width,height= 1, text="Blender",font=self.headerFont,anchor="n", bg=self.canvasbackground, fg =self.canvasforeground).grid(row=12, column=0, columnspan=2, pady=self.separate_label_pad_y)
        self.blenderButton = Button(self.FrameButton,width = self.btn_double_width , height = self.btn_height , text="Activate", font=self.btnFont, command=blendButton, fg=self.btnforeground, bg=self.btnbackground, activebackground=self.btn_common_bg, activeforeground=self.btnforeground)
        self.blenderButton.grid(row=13 ,column=0, columnspan=2,padx=self.pad_x_right,pady=(10,0), sticky="N")
        CreateToolTip(self.blenderButton, self.ttip_dict['blender'])
        self.blend_slider = Scale(self.FrameButton, from_=0, to=1, resolution=0.05, orient=HORIZONTAL, command=self.blendSlices,length = self.scale_length, bg=self.canvasbackground, fg =self.canvasforeground )
        self.blend_slider.grid(row=14,column=0, columnspan = 2, sticky =NW, padx=self.pad_x_left,pady=self.separate_label_pad_y)
        self.blend_slider.set(0.5)
        self.blender = False
        
    def callToolsMenu():
        """ 
            Callback when Tool button clicked.
            
            Constructing Tools Menu.
        """ 
        if self.FrameButton.winfo_exists():
            self.FrameButton.destroy()
        
        self.FrameButton = Frame(self.master, width = self.menuwidth , height = self.option_height, bg=self.canvasbackground)
        self.FrameButton.grid(row = 0, column = 1, rowspan=2, sticky='nw')
        
        createStatusLabel()
        
        Label(self.FrameButton, width = self.btn_single_width,height= 1, text="Resolution",font=self.headerFont,anchor="n", bg=self.canvasbackground, fg =self.canvasforeground).grid(row=0, column=0, columnspan=2, pady=self.frame_pad_y)
        
        Label(self.FrameButton, width = 20,height= 1, text="Transversal [microns]:",font=self.textFont,anchor="nw", bg=self.canvasbackground, fg =self.canvasforeground).grid(row=1,column=0, columnspan=2, padx=self.pad_x_left,pady=(5,0), sticky="W")  
        self.entryTransversal = Entry(self.FrameButton, textvariable=self.tv_res,font=self.textFont)
        self.entryTransversal.grid(row=1,column=0, columnspan=2,padx = self.entry_pad_x,pady=(5,0), sticky="NW")
        CreateToolTip(self.entryTransversal, self.ttip_dict['transversal'])
        
        Label(self.FrameButton, width = 20,height= 1, text="Axial [microns]:",font=self.textFont,anchor="nw", bg=self.canvasbackground, fg =self.canvasforeground).grid(row=2,column=0, columnspan=2, padx=self.pad_x_left,pady=(5,0), sticky="W")  
        self.entryAxial = Entry(self.FrameButton, textvariable=self.axial_res,font=self.textFont)
        self.entryAxial.grid(row=2,column=0, columnspan=2,padx =  self.entry_pad_x,pady=(5,0), sticky="NW")
        CreateToolTip(self.entryAxial, self.ttip_dict['axial'])
        
        self.applyResolution = Button(self.FrameButton,width = self.btn_double_width , height = self.btn_height , text="Activate",font=self.btnFont, command=lambda: background(updateMetrics), fg=self.btnforeground, bg=self.btnbackground, activebackground=self.btn_common_bg, activeforeground=self.btnforeground)
        self.applyResolution.grid(row=3 ,column=0, columnspan=2,padx=self.pad_x_right,pady=(10,0), sticky="N")
        CreateToolTip(self.applyResolution, self.ttip_dict['activateRes'])
 
        Label(self.FrameButton, width = self.btn_single_width,height= 1, text="Metrics",font=self.headerFont,anchor="n", bg=self.canvasbackground, fg =self.canvasforeground).grid(row=4, column=0, columnspan=2, pady=self.separate_label_pad_y)
        
        Label(self.FrameButton, width = 20,height= 1, text="Area of Drusen [mm\u00b2]:",font=self.textFont,anchor="nw", bg=self.canvasbackground, fg =self.canvasforeground).grid(row=5, column=0, pady=self.separate_label_pad_y)
        self.labelDrusenArea = Label(self.FrameButton, width = 20,height= 1, text="0",font=self.textFont,anchor="nw", bg=self.canvasbackground, fg =self.canvasforeground)
        self.labelDrusenArea.grid(row=5, column=1, pady=self.separate_label_pad_y)
        
        Label(self.FrameButton, width = 20,height= 1, text="Area of Drusen [%]:",font=self.textFont,anchor="nw", bg=self.canvasbackground, fg =self.canvasforeground).grid(row=6, column=0, pady=(10,0))
        self.labelDrusenAreaPercent = Label(self.FrameButton, width = 20,height= 1, text="0",font=self.textFont,anchor="nw", bg=self.canvasbackground, fg =self.canvasforeground)
        self.labelDrusenAreaPercent.grid(row=6, column=1, pady=(10,0))

        Label(self.FrameButton, width = 20,height= 1, text="Volume of Drusen [mm\u00b3]:",font=self.textFont,anchor="nw", bg=self.canvasbackground, fg =self.canvasforeground).grid(row=7, column=0, pady=(10,0))
        self.labelDrusenVolume = Label(self.FrameButton, width = 20,height= 1, text="0",font=self.textFont,anchor="nw", bg=self.canvasbackground, fg =self.canvasforeground)
        self.labelDrusenVolume.grid(row=7, column=1, pady=(10,0))
        
        Label(self.FrameButton, width = 20,height= 1, text="Volume of Drusen [%]:",font=self.textFont,anchor="nw", bg=self.canvasbackground, fg =self.canvasforeground).grid(row=8, column=0, pady=(10,0))
        self.labelDrusenVolumePercent = Label(self.FrameButton, width = 20,height= 1, text="0",font=self.textFont,anchor="nw", bg=self.canvasbackground, fg =self.canvasforeground)
        self.labelDrusenVolumePercent.grid(row=8, column=1, pady=(10,0))
        updateMetrics()

    def callDataOptions(event):
        """ 
            Callback when Data button clicked.
            
            Optional
            ----------
            event: event
                unused event
        """ 
        if self.corr_active is False:
            self.canvasData.itemconfig(self.imgData, image = self.imageDataIconActive)
            self.canvasImageOptions.itemconfig(self.imgCon, image = self.imageConIconInactive)
            self.canvasSegmentation.itemconfig(self.imgSeg, image = self.imageSegIconInactive)
            self.canvasVisualization.itemconfig(self.imgVis, image = self.imageVisIconInactive)
            self.canvasTools.itemconfig(self.imgVis, image = self.imageToolsIconInactive)
            callDataMenu()
    
    def callImageOptions(event):
        """ 
            Callback when Image enhancement button clicked.
            
            Optional
            ----------
            event: event
                unused event
        """ 
        if self.corr_active is False:
            self.canvasData.itemconfig(self.imgData, image = self.imageDataIconInactive)
            self.canvasImageOptions.itemconfig(self.imgCon, image = self.imageConIconActive)
            self.canvasSegmentation.itemconfig(self.imgSeg, image = self.imageSegIconInactive)
            self.canvasVisualization.itemconfig(self.imgVis, image = self.imageVisIconInactive)
            self.canvasTools.itemconfig(self.imgVis, image = self.imageToolsIconInactive)
            callImageMenu()
        
    def callSegOptions(event):
        """ 
            Callback when Segmentation button clicked.

            Optional
            ----------
            event: event
                unused event
        """ 
        if self.corr_active is False:
            self.canvasData.itemconfig(self.imgData, image = self.imageDataIconInactive)
            self.canvasImageOptions.itemconfig(self.imgCon, image = self.imageConIconInactive)
            self.canvasSegmentation.itemconfig(self.imgSeg, image = self.imageSegIconActive)
            self.canvasVisualization.itemconfig(self.imgVis, image = self.imageVisIconInactive)
            self.canvasTools.itemconfig(self.imgVis, image = self.imageToolsIconInactive)
            callSegmentationMenu()
        
    def callVisOptions(event):
        """ 
            Callback when Visualization button clicked.
            
            Optional
            ----------
            event: event
                unused event
        """ 
        if self.corr_active is False:
            self.canvasData.itemconfig(self.imgData, image = self.imageDataIconInactive)
            self.canvasImageOptions.itemconfig(self.imgCon, image = self.imageConIconInactive)
            self.canvasSegmentation.itemconfig(self.imgSeg, image = self.imageSegIconInactive)
            self.canvasVisualization.itemconfig(self.imgVis, image = self.imageVisIconActive)
            self.canvasTools.itemconfig(self.imgVis, image = self.imageToolsIconInactive)
            callVisualizationMenu()
            
    def callToolsOptions(event):
        """ 
            Callback when tool button clicked.

            Optional
            ----------
            event: event
                unused event
        """ 
        if self.corr_active is False:
            self.canvasData.itemconfig(self.imgData, image = self.imageDataIconInactive)
            self.canvasImageOptions.itemconfig(self.imgCon, image = self.imageConIconInactive)
            self.canvasSegmentation.itemconfig(self.imgSeg, image = self.imageSegIconInactive)
            self.canvasVisualization.itemconfig(self.imgVis, image = self.imageVisIconInactive)
            self.canvasTools.itemconfig(self.imgVis, image = self.imageToolsIconActive)
            callToolsMenu()
            
    def im_loader(im_string, im_w, im_h):
        """ 
            Load image from string and if small window is true, resize it.
            
            Parameters
            ----------
            im_string: string
                path to image
            im_w: scalar
                image width
            im_h: scalar
                image height
        """ 
        if self.SMALLWINDOW is True:
            return Image.open(os.path.dirname(os.path.abspath(sys.argv[0]))+"\\"+im_string).resize((im_w, im_h), Image.ANTIALIAS)
        else:
            return Image.open(os.path.dirname(os.path.abspath(sys.argv[0]))+"\\"+im_string) 
        
    '''
        Creating basic menu and loading icons to canvas
    '''
    #Basic GUI layout options
    if self.SMALLWINDOW is True:
        heightY = 630 #complete menu height
        im_w = 48 #icon width
        im_h = 125 #icon height
        self.menuwidth = 150 #menu width
        self.option_height = 600 #menu height
        self.btn_double_width = 38 # double column button width
        self.btn_single_width = 14 # single column button width
        self.btn_height = 1 # button height
        self.pad_x_left = (30,0) # single column button padding horizontal (left col)
        self.pad_x_right = (30,53)# single column button padding horizontal (right col)
        self.scale_length = 268 #slider lengths
        self.entry_pad_x=(178,0) #entry widget padding
        self.frame_pad_y = (75,0) # frame vertical padding
        self.separate_label_pad_y=(10,0) # label vertical padding
    else:
        heightY = 980 #menu height
        im_w = 75 #icon width
        im_h = 200 #icon height
        self.menuwidth = 250 #menu width
        self.option_height = 800  #menu height
        self.btn_double_width = 38 # double column button width
        self.btn_single_width = 15 # single column button width
        self.btn_height = 2 # button height
        self.pad_x_left = (40,0) # single column button padding horizontal (left col)
        self.pad_x_right = (40,53) # single column button padding horizontal (right col)
        self.scale_length = 344 #slider lengths
        self.entry_pad_x=(207,0) #entry widget padding
        self.frame_pad_y = (100,0) # frame vertical padding
        self.separate_label_pad_y=(25,0) # label vertical padding

    self.FrameOptions = Frame(self.master, width = im_w , height = heightY, bg=self.canvasbackground)
    self.FrameOptions.grid(row = 0, column = 0, rowspan = 2, sticky='nw')
    
    #Data canvas
    self.canvasData= Canvas(self.FrameOptions, width = im_w, height = im_h, bg=self.canvasbackground,highlightthickness = 0)
    self.canvasData.grid(row=0,column=0, sticky=NW)

    self.canvasData.bind("<Button-1>", callDataOptions)
    img = im_loader('icons//Data_inactive.png',im_w,im_h)
    self.imageDataIconInactive = ImageTk.PhotoImage(img)
    img = im_loader('icons//Data_active.png',im_w,im_h)
    self.imageDataIconActive = ImageTk.PhotoImage(img)
    self.imgData = self.canvasData.create_image(0, 0, image=self.imageDataIconInactive,anchor=NW)
    #Image options canvas
    self.canvasImageOptions = Canvas(self.FrameOptions, width = im_w, height = im_h, bg=self.canvasbackground,highlightthickness = 0)
    self.canvasImageOptions.grid(row=1,column=0)
    self.canvasImageOptions.bind("<Button-1>", callImageOptions)
    img = im_loader('icons//Con_inactive.png',im_w,im_h)
    self.imageConIconInactive = ImageTk.PhotoImage(img)
    img = im_loader('icons//Con_active.png',im_w,im_h)
    self.imageConIconActive = ImageTk.PhotoImage(img)
    self.imgCon = self.canvasImageOptions.create_image(0, 0, image=self.imageConIconInactive,anchor=NW)
    #Segmentation canvas
    self.canvasSegmentation = Canvas(self.FrameOptions, width = im_w, height = im_h, bg=self.canvasbackground,highlightthickness = 0)
    self.canvasSegmentation.grid(row=2,column=0)
    self.canvasSegmentation.bind("<Button-1>", callSegOptions)
    img = im_loader('icons//Seg_inactive.png',im_w,im_h)
    self.imageSegIconInactive = ImageTk.PhotoImage(img)
    img = im_loader('icons//Seg_active.png',im_w,im_h)
    self.imageSegIconActive = ImageTk.PhotoImage(img)
    self.imgSeg = self.canvasSegmentation.create_image(0, 0, image=self.imageSegIconInactive,anchor=NW)
    #Visualization canvas
    self.canvasVisualization = Canvas(self.FrameOptions, width = im_w, height = im_h, bg=self.canvasbackground,highlightthickness = 0)
    self.canvasVisualization.grid(row=3,column=0)
    self.canvasVisualization.bind("<Button-1>", callVisOptions)
    img = im_loader('icons//Vis_inactive.png',im_w,im_h)
    self.imageVisIconInactive = ImageTk.PhotoImage(img)
    img = im_loader('icons//Vis_active.png',im_w,im_h)
    self.imageVisIconActive = ImageTk.PhotoImage(img)
    self.imgVis = self.canvasVisualization.create_image(0, 0, image=self.imageVisIconInactive,anchor=NW)
    #Tools canvas
    self.canvasTools = Canvas(self.FrameOptions, width = im_w, height = im_h, bg=self.canvasbackground,highlightthickness = 0)
    self.canvasTools.grid(row=4,column=0)
    self.canvasTools.bind("<Button-1>", callToolsOptions)
    img = im_loader('icons//Tools_inactive.png',im_w,im_h)
    self.imageToolsIconInactive = ImageTk.PhotoImage(img)
    img = im_loader('icons//Tools_active.png',im_w,im_h)
    self.imageToolsIconActive = ImageTk.PhotoImage(img)
    self.imgTools = self.canvasTools.create_image(0, 0, image=self.imageToolsIconInactive,anchor=NW)

    self.FrameButton = Frame(self.master, width = self.menuwidth , height = self.option_height, bg=self.canvasbackground)
    self.FrameButton.grid(row = 0, column = 1, rowspan=2, sticky='nw')
    
    #call segmentation when coming from correction, else data
    if(self.fromCorrMenu):
        callSegOptions(None)
    else:
        callDataOptions(None)