""" 
CanvasHeatmap: Module to create heatmap canvas
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
import PIL.Image, PIL.ImageTk
from tkinter import *
from matplotlib import cm
from Visualization import Heatmap
from Tooltips.TipHandler import CreateToolTip

def createCanvasHeatmap(self):
    """ 
        Creation of heatmap and event handlers.
        
        The colormap is configured by Image Menu-> Colormap Dropdown.
        
        Events:
        Right-click on canvas: show intensity at position (x,y)
        Mouse-wheel: Threshold slider up/down
        
        Parameters:
        ----------
        self: object
            Framework
    """
    def mousewheelHeatmap(event):
        if event.num == 5 or event.delta == -120:
            self.heat_slider_min.set(self.heat_slider_min.get()-1)
        if event.num == 4 or event.delta == 120:
            self.heat_slider_min.set(self.heat_slider_min.get()+1)
    
    def rightclickHeatmap(event):
        if self.SMALLWINDOW is True:
            self.LabelIntensityHM.configure(text = str(int(self.heatmapThreshed[int(event.y*2/self.scale_img),int(event.x*2/self.scale_img)])))
        else:
            self.LabelIntensityHM.configure(text = str(int(self.heatmapThreshed[int(event.y/self.scale_img),int(event.x/self.scale_img)])))
    
    self.heatmapThreshed = 0
    self._job_HM = None
    #calculate heatmap and set all heatmaps
    self.heatmap_rpe_bruchs = Heatmap.calculateHeatmap(self.segmentation)
    self.heatmap = self.heatmap_rpe_bruchs 
    self.heatmapThreshed = self.heatmap_rpe_bruchs 
    
    #setup frame
    self.heightHM, self.widthHM = self.heatmap.shape
    if(self.SMALLWINDOW is True):
        self.FrameHeatmap = Frame(self.master, width = 1, height = 1, relief=SUNKEN, bg=self.canvasbackground)
    else:
        self.FrameHeatmap = Frame(self.master, width = 1, height = 1, bg=self.canvasbackground,relief=SUNKEN)
    self.FrameHeatmap.grid(row= 0, column = 3,sticky ='nw')
    
    #Plane and Label for heatmap
    if(self.SMALLWINDOW is True):
        self.LabelIntensityHM = Label(self.FrameHeatmap, text = 'Right click on pixel',width =31, font=self.clrbarFont, bg=self.canvasbackground,fg =self.canvasforeground)
    else:
        self.LabelIntensityHM = Label(self.FrameHeatmap, text = 'Right click on pixel',width =62, font=self.clrbarFont, bg=self.canvasbackground,fg =self.canvasforeground)
 
    self.LabelIntensityHM.grid(row = 0, column = 0)
    CreateToolTip(self.LabelIntensityHM, self.ttip_dict['heatmapIntensity'])
    if(self.SMALLWINDOW is True):
        self.canvas_heat = Canvas(self.FrameHeatmap, width = self.widthHM//2, height = self.heightHM//2, bg=self.canvasbackground,highlightthickness=2,highlightbackground=self.hlghtbg)
    else:
        self.canvas_heat = Canvas(self.FrameHeatmap, width = self.widthHM, height = self.heightHM, bg=self.canvasbackground,highlightthickness=2,highlightbackground=self.hlghtbg)
    self.canvas_heat.grid(row=1,column=0, rowspan=3)
    
    #mousebindings
    self.canvas_heat.bind("<Button-3>", rightclickHeatmap)
    self.canvas_heat.bind("<MouseWheel>", mousewheelHeatmap)

    #scrollbars
    self.sbarV_HM = Scrollbar(self.FrameHeatmap, orient=VERTICAL,command=self.getTopYMotion)
    self.sbarH_HM = Scrollbar(self.FrameHeatmap, orient=HORIZONTAL,command=self.getTopXMotion)
    self.canvas_heat.config(yscrollcommand=self.sbarV_HM.set)
    self.canvas_heat.config(xscrollcommand=self.sbarH_HM.set)
    self.sbarV_HM.grid(row=1, column = 1, sticky=N+S, rowspan=3)
    self.sbarH_HM.grid(row=4, column = 0, sticky=W+E)
    self.canvas_heat.configure(scrollregion = self.canvas_heat.bbox("all"))
    # Colorbar
    self.heatmap_max = np.maximum(10,np.max(self.heatmap))
    
    self.label_top_HM = Label(self.FrameHeatmap, width = 20,height= 1, text=' '+str(int(self.heatmap_max)),font=self.clrbarFont,anchor="nw", bg=self.canvasbackground, fg =self.canvasforeground)
    self.label_top_HM.grid(row=1,column=2)
    if(self.SMALLWINDOW == True):
        self.canvas_cbar_HM = Canvas(self.FrameHeatmap, width = 18, height = self.heightHM//2 - 50, bg=self.canvasbackground)
    else:
        self.canvas_cbar_HM = Canvas(self.FrameHeatmap, width=18, height = self.heightHM - 50, bg=self.canvasbackground)
    self.canvas_cbar_HM.grid(row=2,column=2,padx=(3, 0),pady = (0,0),sticky="nw")
    self.label_bottom_HM = Label(self.FrameHeatmap, width = 20, height=  1, text="  0",font=self.clrbarFont,anchor="nw", bg=self.canvasbackground, fg =self.canvasforeground)
    self.label_bottom_HM.grid(row=3,column=2)
    
    self.cmap_own= cm.get_cmap(self.COLORMAP, self.heatmap_max)
    
    if(self.SMALLWINDOW == True):
        
        cmap_std= cm.get_cmap(self.COLORMAP, 256)
            
        self.cbar_arr = np.array([range(0,self.heightHM//2)[::-1] for i in range(20)]).T.astype('uint8')
        self.heatmapThreshedPil = PIL.Image.fromarray(self.cmap_own(self.heatmap[::2,::2], bytes=True))
    else:
        cmap_std= cm.get_cmap(self.COLORMAP, 256)

        self.cbar_arr = np.zeros((self.heightHM, 20)).astype('uint8')
        
        for i in range (256):
            self.cbar_arr[i*2:i*2+2,:] = i
        self.cbar_arr = self.cbar_arr[::-1]
        self.heatmapThreshedPil = PIL.Image.fromarray(self.cmap_own(self.heatmap, bytes=True))
    
    #update canvas with images
    self.photo_heat = PIL.ImageTk.PhotoImage(image = self.heatmapThreshedPil)
    self.photo_cbar = PIL.ImageTk.PhotoImage(image = PIL.Image.fromarray(cmap_std(self.cbar_arr, bytes=True)))
    
    self.image_on_canvas_heat = self.canvas_heat.create_image(0, 0, image=self.photo_heat, anchor=NW)
    self.image_on_canvas_cbar = self.canvas_cbar_HM.create_image(0, 0, image=self.photo_cbar, anchor=NW)
    
    if(self.SMALLWINDOW == True):
        self.heat_slider_min = Scale(self.FrameHeatmap, from_=0, to=self.heatmap_max, orient=HORIZONTAL, length = self.widthHM//2, command=self.updateHeatSlice,bg=self.canvasbackground, fg =self.canvasforeground)
    else:
        self.heat_slider_min = Scale(self.FrameHeatmap, from_=0, to=self.heatmap_max, orient=HORIZONTAL, length = self.widthHM, command=self.updateHeatSlice,bg=self.canvasbackground, fg =self.canvasforeground)
    self.heat_slider_min.grid(row=5,column=0)
    self.heat_slider_min.set(0)
    CreateToolTip(self.heat_slider_min, self.ttip_dict['heatmapThresh'])