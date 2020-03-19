""" 
CanvasYZ: Module to create yz Canvas
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

def createCanvasYZ(self):
    """ 
        Creation of YZ plane and configuring event handlers.
        
        Events:
        Mouse-wheel: Scroll through slices
        Left-click: Set all planes to specific point and show crosshair (2s)
        Right-click: Show intensity at specific location
    
        Parameters:
        ----------
        self: object
            Framework
    """
    def mousewheelYZ(event):
        """ 
            Mouse wheel yz - change slice
        """
        if event.num == 5 or event.delta == -120:
            self.vol_sliderYZ.set(self.vol_sliderYZ.get()-1)
        if event.num == 4 or event.delta == 120:
            self.vol_sliderYZ.set(self.vol_sliderYZ.get()+1)
            
    def leftclickYZ(event):
        """ 
            Left click yz - creates cross hair on planes for 2000 ms
        """
        z = self.vol_sliderYZ.get()
        
        canvas_event = event.widget
        c_x = int(canvas_event.canvasx(event.x))
        c_y = int(canvas_event.canvasy(event.y))

        if(self.SMALLWINDOW is True):  
            self.cross_hz_xy = self.canvasXY.create_line(self.scale_img*z//2, 0, self.scale_img*z//2, self.scale_img*self.widthXY//2, fill=self.crosshair_color, width=1)
            self.cross_vt_xy = self.canvasXY.create_line(0, c_x, self.scale_img*self.heightXY//2, c_x, fill=self.crosshair_color, width=1)
            self.cross_hz_xz = self.canvasXZ.create_line(self.scale_img*z//2, 0, self.scale_img*z//2, self.scale_img*self.heightXZ//2, fill=self.crosshair_color, width=1)
            self.cross_vt_xz = self.canvasXZ.create_line(0, c_y, self.scale_img*self.widthXZ//2, c_y, fill=self.crosshair_color, width=1)
            self.cross_hz_yz = self.canvasYZ.create_line(c_x, 0, c_x, self.scale_img*self.heightYZ//2, fill=self.crosshair_color, width=1)
            self.cross_vt_yz = self.canvasYZ.create_line(0, c_y, self.scale_img*self.widthYZ//2, c_y, fill=self.crosshair_color, width=1)
            c_x = c_x * 2
            c_y = c_y * 2
        else:
            self.cross_hz_xy = self.canvasXY.create_line(self.scale_img*z, 0, self.scale_img*z, self.scale_img*self.widthXY, fill=self.crosshair_color, width=1)
            self.cross_vt_xy = self.canvasXY.create_line(0, c_x, self.scale_img*self.heightXY, c_x, fill=self.crosshair_color, width=1)
            self.cross_hz_xz = self.canvasXZ.create_line(self.scale_img*z, 0, self.scale_img*z, self.scale_img*self.heightXZ, fill=self.crosshair_color, width=1)
            self.cross_vt_xz = self.canvasXZ.create_line(0, c_y, self.scale_img*self.widthXZ, c_y, fill=self.crosshair_color, width=1)
            self.cross_hz_yz = self.canvasYZ.create_line(c_x, 0, c_x, self.scale_img*self.heightYZ, fill=self.crosshair_color, width=1)
            self.cross_vt_yz = self.canvasYZ.create_line(0, c_y, self.scale_img*self.widthYZ, c_y, fill=self.crosshair_color, width=1)

        self.canvasXZ.after(2000, self.canvasXZ.delete, self.cross_hz_xz)
        self.canvasXZ.after(2000, self.canvasXZ.delete, self.cross_vt_xz)
        self.canvasXY.after(2000, self.canvasXY.delete, self.cross_hz_xy)
        self.canvasXY.after(2000, self.canvasXY.delete, self.cross_vt_xy)
        self.canvasYZ.after(2000, self.canvasYZ.delete, self.cross_hz_yz)
        self.canvasYZ.after(2000, self.canvasYZ.delete, self.cross_vt_yz)
        self.vol_sliderXY.set(c_y//self.scale_img)
        self.vol_sliderXZ.set(c_x//self.scale_img)
        
    def rightclickYZ(event):
        """ 
            Right click YZ - gets intensity at pixel or error
        """
        if self.SMALLWINDOW == True:
            try:
                self.LabelDistanceXY.configure(text = str(int(self.volume_original[int(event.x*2/self.scale_img),int(event.y*2/self.scale_img),self.vol_sliderYZ.get()])))
            except:
                self.LabelDistanceXY.configure(text = "Error")
        else:
            try:
                self.LabelDistanceXY.configure(text = str(int(self.volume_original[int(event.x/self.scale_img),int(event.y/self.scale_img),self.vol_sliderYZ.get()])))
            except:
                self.LabelDistanceXY.configure(text = "Error")
    ''' 
        Initial creation process
    '''            
    self._job_YZ = None
    
    self.sliceYZ = np.zeros((self.volume.shape[0],self.volume.shape[2], 3)).astype('uint8')
    
    self.heightYZ = self.volume.shape[0]
    self.widthYZ = self.volume.shape[2]
    #basic frame
    if(self.SMALLWINDOW == True):
        self.FrameVolumeYZ = Frame(self.master, width = self.widthYZ//2 , height = self.heightYZ//2 + 20, relief=SUNKEN)
    else:
        self.FrameVolumeYZ = Frame(self.master, width = self.widthYZ , height = 350 + 20, relief=SUNKEN)
    self.FrameVolumeYZ.grid(row = 1, column = 3,sticky ='nw')
    #canvas
    if(self.SMALLWINDOW == True):    
        self.canvasYZ = Canvas(self.FrameVolumeYZ, width = self.widthYZ//2, height = self.heightYZ//2,highlightthickness=2,highlightbackground=self.hlghtbg)
    else:
        self.canvasYZ = Canvas(self.FrameVolumeYZ, width = self.widthYZ, height = 350,highlightthickness=2,highlightbackground=self.hlghtbg)
    self.canvasYZ.grid(row=0,column=0)
    #mouse bindings
    self.canvasYZ.bind("<Button-1>", leftclickYZ)
    self.canvasYZ.bind("<Button-3>", rightclickYZ)
    self.canvasYZ.bind("<MouseWheel>", mousewheelYZ)
    #scrollbar
    self.sbarV_YZ = Scrollbar(self.FrameVolumeYZ, orient=VERTICAL,command=self.getBottomYMotion)
    self.sbarH_YZ = Scrollbar(self.FrameVolumeYZ, orient=HORIZONTAL,command=self.getBottomXMotion)
    self.canvasYZ.config(yscrollcommand=self.sbarV_YZ.set)
    self.canvasYZ.config(xscrollcommand=self.sbarH_YZ.set)
    self.sbarV_YZ.grid(row=0, column = 1, sticky=N+S)
    self.sbarH_YZ.grid(row=1, column = 0, sticky=W+E)
    self.canvasYZ.configure(scrollregion = self.canvasYZ.bbox("all"))
    
    #update planes
    self.volumeYZ = np.swapaxes(self.volume, axis1 = 0, axis2 = 1)
    self.volumeYZ = np.swapaxes(self.volumeYZ, axis1 = 0, axis2 = 2)
    
    if(self.SMALLWINDOW == True): 
        self.slice_pilYZ = PIL.Image.fromarray(self.sliceYZ[::2,::2])
    else:
        self.slice_pilYZ = PIL.Image.fromarray(self.sliceYZ)
    self.photoYZ = PIL.ImageTk.PhotoImage(image = self.slice_pilYZ)
    
    self.image_on_canvas_YZ = self.canvasYZ.create_image(0, 0, image=self.photoYZ, anchor=NW)
    #scale
    if(self.SMALLWINDOW == True):
        self.vol_sliderYZ = Scale(self.FrameVolumeYZ, from_=0, to=self.volumeYZ.shape[0]-1, length = self.widthYZ//2,orient=HORIZONTAL, command=self.updateVolumeYZSlice, bg=self.canvasbackground, fg =self.canvasforeground)
    else:
        self.vol_sliderYZ = Scale(self.FrameVolumeYZ, from_=0, to=self.volumeYZ.shape[0]-1, length = self.widthYZ,orient=HORIZONTAL, command=self.updateVolumeYZSlice, bg=self.canvasbackground, fg =self.canvasforeground)
    self.vol_sliderYZ.grid(row=2,column=0)
    self.vol_sliderYZ.set(self.volumeYZ.shape[0]//2)