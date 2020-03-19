""" 
CanvasXZ: Module to create XZ Canvas
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

def createCanvasXZ(self):
    """ 
        Creation of XZ plane and configurung event handlers.
        
        Events:
        Mouse-wheel: Scroll through slices
        Left-click: Set all planes to specific point and show crosshair (2s)
        Right-click: Show intensity at specific location
    
        Parameters:
        ----------
        self: object
            Framework
    """
    def mousewheelXZ(event):
        """ 
            Mouse wheel XZ - change slice
        """
        if event.num == 5 or event.delta == -120:
            self.vol_sliderXZ.set(self.vol_sliderXZ.get()-1)
        if event.num == 4 or event.delta == 120:
            self.vol_sliderXZ.set(self.vol_sliderXZ.get()+1)
    
    def leftclickXZ(event):
        """ 
            Left click XZ - create crosshair for 2000 ms
        """
        z = self.vol_sliderXZ.get()
        
        canvas_event = event.widget
        c_x = int(canvas_event.canvasx(event.x))
        c_y = int(canvas_event.canvasy(event.y))
        
        if(self.SMALLWINDOW is True):
            self.cross_hz_xy = self.canvasXY.create_line(c_x, 0, c_x, self.scale_img*self.heightXY//2, fill=self.crosshair_color, width=1)
            self.cross_vt_xy = self.canvasXY.create_line(0, self.scale_img*z//2, self.scale_img*self.widthXY//2, z//2, fill=self.crosshair_color, width=1)
            self.cross_hz_xz = self.canvasXZ.create_line(c_x, 0, c_x, self.scale_img*self.heightXZ//2, fill=self.crosshair_color, width=1)
            self.cross_vt_xz = self.canvasXZ.create_line(0, c_y, self.scale_img*self.widthXZ//2, c_y, fill=self.crosshair_color, width=1)
            self.cross_hz_yz = self.canvasYZ.create_line(self.scale_img*z//2,0, self.scale_img*z//2,self.scale_img*self.widthYZ//2, fill=self.crosshair_color, width=1)
            self.cross_vt_yz = self.canvasYZ.create_line(0, c_y, self.scale_img*self.widthYZ//2, c_y, fill=self.crosshair_color, width=1)
            c_x = c_x * 2
            c_y = c_y * 2
        else:
            self.cross_hz_xy = self.canvasXY.create_line(c_x, 0, c_x, self.scale_img*self.heightXY, fill=self.crosshair_color, width=1)
            self.cross_vt_xy = self.canvasXY.create_line(0, self.scale_img*z, self.scale_img*self.widthXY, self.scale_img*z, fill=self.crosshair_color, width=1)
            self.cross_hz_xz = self.canvasXZ.create_line(c_x, 0, c_x, self.scale_img*self.heightXZ, fill=self.crosshair_color, width=1)
            self.cross_vt_xz = self.canvasXZ.create_line(0, c_y, self.scale_img*self.widthXZ, c_y, fill=self.crosshair_color, width=1)
            self.cross_hz_yz = self.canvasYZ.create_line(self.scale_img*z,0, self.scale_img*z,self.scale_img*self.widthYZ, fill=self.crosshair_color, width=1)
            self.cross_vt_yz = self.canvasYZ.create_line(0, c_y, self.scale_img*self.widthYZ, c_y, fill=self.crosshair_color, width=1)
        
        self.canvasYZ.after(2000, self.canvasYZ.delete, self.cross_hz_yz)
        self.canvasYZ.after(2000, self.canvasYZ.delete, self.cross_vt_yz)
        self.canvasXZ.after(2000, self.canvasXZ.delete, self.cross_hz_xz)
        self.canvasXZ.after(2000, self.canvasXZ.delete, self.cross_vt_xz)
        self.canvasXY.after(2000, self.canvasXY.delete, self.cross_hz_xy)
        self.canvasXY.after(2000, self.canvasXY.delete, self.cross_vt_xy)
        self.vol_sliderXY.set(c_y//self.scale_img)
        self.vol_sliderYZ.set(c_x//self.scale_img)
        
    def rightclickXZ(event):
        """ 
            Right click XZ - show intensity in en face label
        """
        if self.SMALLWINDOW == True:
            try:
                self.LabelDistanceXY.configure(text = str(int(self.volume_original[self.vol_sliderXZ.get(),int(event.y*2/self.scale_img),int(event.x*2/self.scale_img)])))
            except:
                self.LabelDistanceXY.configure(text = "Error")
        else:
            try:
                self.LabelDistanceXY.configure(text = str(int(self.volume_original[self.vol_sliderXZ.get(),int(event.y/self.scale_img),int(event.x/self.scale_img)])))
            except:
                self.LabelDistanceXY.configure(text = "Error")
    
    '''
        basic creation process
    '''
    self._job_XZ = None
    
    self.sliceXZ = np.zeros((self.volume.shape[0],self.volume.shape[2],3)).astype('uint8')
    
    self.heightXZ = self.volume.shape[0]
    self.widthXZ = self.volume.shape[2]
    #basic frame
    if(self.SMALLWINDOW == True):    
        self.FrameVolumeXZ = Frame(self.master, width = self.widthXZ//2 , height = self.heightXZ//2 + 20, relief=SUNKEN)
    else:
        self.FrameVolumeXZ = Frame(self.master, width = self.widthXZ , height = 350 + 20, relief=SUNKEN)
    self.FrameVolumeXZ.grid(row = 1, column = 2, sticky ='nw')
    #canvas
    if(self.SMALLWINDOW == True):
        self.canvasXZ = Canvas(self.FrameVolumeXZ, width = self.widthXZ//2, height = self.heightXZ//2,highlightthickness=2,highlightbackground=self.hlghtbg)
    else:
        self.canvasXZ = Canvas(self.FrameVolumeXZ, width = self.widthXZ, height = 350,highlightthickness=2,highlightbackground=self.hlghtbg)
    self.canvasXZ.grid(row=0,column=0)
    
    #mouse bindings    
    self.canvasXZ.bind("<Button-1>", leftclickXZ)
    self.canvasXZ.bind("<Button-3>", rightclickXZ)
    self.canvasXZ.bind("<MouseWheel>", mousewheelXZ)
    #scroll bars
    self.sbarV_XZ = Scrollbar(self.FrameVolumeXZ, orient=VERTICAL,command=self.getBottomYMotion)
    self.sbarH_XZ = Scrollbar(self.FrameVolumeXZ, orient=HORIZONTAL,command=self.getBottomXMotion)
    self.canvasXZ.config(yscrollcommand=self.sbarV_XZ.set)
    self.canvasXZ.config(xscrollcommand=self.sbarH_XZ.set)
    self.sbarV_XZ.grid(row=0, column = 1, sticky=N+S)
    self.sbarH_XZ.grid(row=1, column = 0, sticky=W+E)
    self.canvasXZ.configure(scrollregion = self.canvasXZ.bbox("all"))
    
    #load images to canvas
    self.volumeXZ = np.swapaxes(self.volume, axis1 = 1, axis2 = 0)
    
    if(self.SMALLWINDOW == True):
        slice_pilXZ = PIL.Image.fromarray(self.sliceXZ[::2,::2])
    else:
        slice_pilXZ = PIL.Image.fromarray(self.sliceXZ)

    self.image_on_canvas_XZ = self.canvasXZ.create_image(0, 0, image=PIL.ImageTk.PhotoImage(image = slice_pilXZ), anchor=NW)
    #scale
    if self.SMALLWINDOW == True:
        self.vol_sliderXZ = Scale(self.FrameVolumeXZ, from_=0, to=self.volumeXZ.shape[0]-1, length = self.widthXZ//2, orient=HORIZONTAL, command=self.updateVolumeXZSlice, bg=self.canvasbackground, fg =self.canvasforeground)
    else:
        self.vol_sliderXZ = Scale(self.FrameVolumeXZ, from_=0, to=self.volumeXZ.shape[0]-1, length = self.widthXZ, orient=HORIZONTAL, command=self.updateVolumeXZSlice, bg=self.canvasbackground, fg =self.canvasforeground)
    self.vol_sliderXZ.grid(row=2,column=0)
    self.vol_sliderXZ.set(self.volumeXZ.shape[0]//2)