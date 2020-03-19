""" 
CanvasXY: Module to create en face Canvas
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
from Tooltips.TipHandler import CreateToolTip

def createCanvasXY(self):
    """ 
        Creation of En face plane and configuring event handlers.
        
        Events:
        Mouse-wheel: Scroll through slices
        Left-click: Set all planes to specific point and show crosshair (2s)
        Right-click: Show intensity at specific location
    
        Parameters:
        ----------
        self: object
            Framework
    """
    def mousewheelXY(event):
        """ 
            Mouse wheel en face - change slice
        """
        if event.num == 5 or event.delta == -120:
            self.vol_sliderXY.set(self.vol_sliderXY.get()-1)
        if event.num == 4 or event.delta == 120:
            self.vol_sliderXY.set(self.vol_sliderXY.get()+1)
        
    def leftclickXY(event):
        """ 
            Left click en face - creates cross hair on planes for 2000 ms
        """
        z = self.vol_sliderXY.get()
        canvas_event = event.widget
        c_x = int(canvas_event.canvasx(event.x))
        c_y = int(canvas_event.canvasy(event.y))
        
        if(self.SMALLWINDOW is True):  
            self.cross_hz_xy = self.canvasXY.create_line(c_x, 0, c_x, self.scale_img*self.heightXY//2, fill=self.crosshair_color, width=1)
            self.cross_vt_xy = self.canvasXY.create_line(0, c_y, self.scale_img*self.widthXY//2, c_y, fill=self.crosshair_color, width=1)
            self.cross_hz_xz = self.canvasXZ.create_line(c_x, 0, c_x, self.scale_img*self.heightXZ//2, fill=self.crosshair_color, width=1)
            self.cross_vt_xz = self.canvasXZ.create_line(0, (int(z*self.scale_img)//2), self.scale_img*self.widthXZ//2, (int(z*self.scale_img)//2), fill=self.crosshair_color, width=1)
            self.cross_hz_yz = self.canvasYZ.create_line(c_y, 0, c_y, self.scale_img*self.heightYZ//2, fill=self.crosshair_color, width=1)
            self.cross_vt_yz = self.canvasYZ.create_line(0, (int(z*self.scale_img)//2), self.scale_img*self.widthYZ//2, (int(z*self.scale_img)//2), fill=self.crosshair_color, width=1)
            c_x = c_x * 2
            c_y = c_y * 2
        else:
            self.cross_hz_xy = self.canvasXY.create_line(c_x, 0, c_x, self.scale_img*self.heightXY, fill=self.crosshair_color, width=1)
            self.cross_vt_xy = self.canvasXY.create_line(0, c_y, self.scale_img*self.widthXY, c_y, fill=self.crosshair_color, width=1)
            
            #z = z - self.volume.shape[0] + 1
            self.cross_hz_xz = self.canvasXZ.create_line(c_x, 0, c_x, self.scale_img*self.heightXZ, fill=self.crosshair_color, width=1)
            self.cross_vt_xz = self.canvasXZ.create_line(0, int(z*self.scale_img), self.scale_img*self.widthXZ,  int(z*self.scale_img), fill=self.crosshair_color, width=1)
            self.cross_hz_yz = self.canvasYZ.create_line(c_y, 0, c_y, self.scale_img*self.heightYZ, fill=self.crosshair_color, width=1)
            self.cross_vt_yz = self.canvasYZ.create_line(0,int(z*self.scale_img), self.scale_img*self.widthYZ,  int(z*self.scale_img), fill=self.crosshair_color, width=1)     
            
        self.canvasYZ.after(2000, self.canvasYZ.delete, self.cross_hz_yz)
        self.canvasYZ.after(2000, self.canvasYZ.delete, self.cross_vt_yz)
        self.canvasXY.after(2000, self.canvasXY.delete, self.cross_hz_xy)
        self.canvasXY.after(2000, self.canvasXY.delete, self.cross_vt_xy)
        self.canvasXZ.after(2000, self.canvasXZ.delete, self.cross_hz_xz)
        self.canvasXZ.after(2000, self.canvasXZ.delete, self.cross_vt_xz)
        self.vol_sliderXZ.set(c_y//self.scale_img)
        self.vol_sliderYZ.set(c_x//self.scale_img)
        
    def rightclickXY(event):
        """ 
            Right click en face - gets intensity at pixel or error
        """
        print(self.scale_img)
        if self.SMALLWINDOW is True:
            try:
                self.LabelDistanceXY.configure(text = str(int(self.volume_original[int(event.y*2/self.scale_img),self.vol_sliderXY.get(),int(event.x*2/self.scale_img)])))
            except:
                self.LabelDistanceXY.configure(text = "Error")
        else:
            try:
                self.LabelDistanceXY.configure(text = str(int(self.volume_original[int(event.y/self.scale_img),self.vol_sliderXY.get(),int(event.x/self.scale_img)])))
            except:
                self.LabelDistanceXY.configure(text = "Error")
    
    ''' 
        Initial creation process
    '''
    self._job_volume_xy = None
    
    self.slicexy = np.zeros((self.volume.shape[1],self.volume.shape[2], 3)).astype('uint8')
    self.heightXY, self.widthXY = self.slicexy.shape[0:2]
    #basic frame
    if(self.SMALLWINDOW is True):
        self.FrameVolumeXY = Frame(self.master, width = self.widthXY//2 , height = self.heightXY//2 + 40, relief=SUNKEN)
    else:
        self.FrameVolumeXY = Frame(self.master, width = self.widthXY , height = self.heightXY + 40, relief=SUNKEN)
    self.FrameVolumeXY.grid(row = 0, column = 2)
    #top label
    if(self.SMALLWINDOW is True):
        self.LabelDistanceXY = Label(self.FrameVolumeXY, text = 'Right click on pixel', width =35, font=self.clrbarFont,anchor='n', bg=self.canvasbackground, fg =self.canvasforeground)
    else:
        self.LabelDistanceXY = Label(self.FrameVolumeXY, text = 'Right click on pixel', width =62, font=self.clrbarFont,anchor='n', bg=self.canvasbackground, fg =self.canvasforeground)
    self.LabelDistanceXY.grid(row = 0, column = 0)
    CreateToolTip(self.LabelDistanceXY, self.ttip_dict['enfaceIntensity'])
    #en face canvas
    if(self.SMALLWINDOW == True):
        self.canvasXY = Canvas(self.FrameVolumeXY, width = self.widthXY//2, height = self.heightXY//2,highlightthickness=2,highlightbackground=self.hlghtbg)
    else:
        self.canvasXY = Canvas(self.FrameVolumeXY, width = self.widthXY, height = self.heightXY,highlightthickness=2,highlightbackground=self.hlghtbg)
    self.canvasXY.grid(row=1,column=0)
    #mouse bindings
    self.canvasXY.bind("<Button-1>",   leftclickXY)
    self.canvasXY.bind("<Button-3>",   rightclickXY)
    self.canvasXY.bind("<MouseWheel>", mousewheelXY)

    #load images to canvas
    if(self.SMALLWINDOW == True):
        slice_pil_xy = PIL.Image.fromarray(self.slicexy[::2,::2])
    else:
        slice_pil_xy = PIL.Image.fromarray(self.slicexy)
    
    self.image_on_canvas_xy = self.canvasXY.create_image(0, 0, image=PIL.ImageTk.PhotoImage(image = slice_pil_xy), anchor=NW)
    
    self.sbarV_XY = Scrollbar(self.FrameVolumeXY, orient=VERTICAL,command=self.getTopYMotion)
    self.sbarH_XY = Scrollbar(self.FrameVolumeXY, orient=HORIZONTAL,command=self.getTopXMotion)
    self.canvasXY.config(yscrollcommand=self.sbarV_XY.set)
    self.canvasXY.config(xscrollcommand=self.sbarH_XY.set)
    self.sbarV_XY.grid(row=1, column = 1, sticky=N+S)
    self.sbarH_XY.grid(row=2, column = 0, sticky=W+E)
    self.canvasXY.configure(scrollregion = self.canvasXY.bbox("all"))
    
    #scale to scroll through volume
    if(self.SMALLWINDOW == True):
        self.vol_sliderXY = Scale(self.FrameVolumeXY, from_=0, to=self.volume.shape[0]-1, length = self.widthXY//2, orient=HORIZONTAL, command=self.updateVolumeSliceXY, bg=self.canvasbackground, fg =self.canvasforeground)
    else:
        self.vol_sliderXY = Scale(self.FrameVolumeXY, from_=0, to=self.volume.shape[0]-1, length = self.widthXY, orient=HORIZONTAL, command=self.updateVolumeSliceXY, bg=self.canvasbackground, fg =self.canvasforeground)
    self.vol_sliderXY.grid(row=3,column=0)
    self.vol_sliderXY.set(self.volume.shape[0]//2)