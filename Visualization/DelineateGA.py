""" 
DelineateGA: Module to load GA lesion delineation
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
import cv2

def pack(self):
    """
        Method to delineate GA and calculate covered lesion area.
        
        To automatically delineate GA by a rectangular input:
        1) Left-Click and hold at starting point for rectangle
        2) Release rectangle at end-point
        3) The algorithm automatically calculates the GA area
        
        Note: Works only with Open CV as it uses the Grab-Cut algorithm.
    """
    
    def startRect(event):
        """
            Start rectanglular selection
        """
        
        self.move = True
        self.rectx0 = self.canvas_heat.canvasx(event.x)
        self.recty0 = self.canvas_heat.canvasy(event.y)
        self.rect = self.canvas_heat.create_rectangle(self.rectx0, self.recty0, self.rectx0, self.recty0, outline='#006400',width=3, tags='rectangleGA')
        self.rectid = self.canvas_heat.find_closest(self.rectx0, self.recty0, halo=2)
    
    def movingRect(event):
        """
            Move rectanglular selection
        """
        
        if self.move: 
            self.rectx1 = self.canvas_heat.canvasx(event.x)
            self.recty1 = self.canvas_heat.canvasy(event.y)
            self.canvas_heat.coords(self.rectid, self.rectx0, self.recty0, self.rectx1, self.recty1)
    
    def stopRect(event):
        """
            Stop rectanglular selection and calculate GA lesion
        """
        
        self.move = False

        self.rectx1 = self.canvas_heat.canvasx(event.x)
        self.recty1 = self.canvas_heat.canvasy(event.y) 
        self.canvas_heat.coords(self.rectid, self.rectx0, self.recty0, self.rectx1, self.recty1)
        
        box = (int(self.rectx0), int(self.recty0), int(self.rectx1)-int(self.rectx0), int(self.recty1)-int(self.recty0))
    
        markGaByRectangle(box)
            
    def createVariables(parent):
        """
            Create variables
        """
        self.parent = parent
        self.rectx0 = 0
        self.recty0 = 0
        self.rectx1 = 0
        self.recty1 = 0
        self.rectid = None
        self.move = False
    
    def createCanvasBinding():
        """
            Create mouse-bindings
        """
        self.canvas_heat.bind( "<Button-1>", startRect)
        self.canvas_heat.bind( "<ButtonRelease-1>", stopRect)
        self.canvas_heat.bind( "<Motion>", movingRect)

    def markGaByRectangle(rect):
        """
            Calculate lesion area by using grab-cut
        """
        
        #clear selection
        clearMarkingData(self)
        #rgb image
        img = cv2.cvtColor(self.heatmapThreshed, cv2.COLOR_GRAY2RGB)
        mask = np.zeros(img.shape[:2],np.uint8)
        #generate rectangle if inverted selection
        if(rect[2] < 0):
            r0 = rect[0]+rect[2]
            r2 = -rect[2]
            rect = (r0,rect[1],r2,rect[3])
            
        if(rect[3] < 0):
            r1 = rect[1]+rect[3]
            r3 = -rect[3]
            rect = (rect[0],r1,rect[2],r3)
        #run algorithm
        try:
            #grab-cut
            cv2.grabCut(img, mask, rect, np.zeros((1,65),np.float64), np.zeros((1,65),np.float64), 5, cv2.GC_INIT_WITH_RECT)
            mask2 = np.where((mask==2)|(mask==0),0,1).astype('uint8')
            _, contours, _ = cv2.findContours(mask2,cv2.RETR_LIST,cv2.CHAIN_APPROX_SIMPLE)
            #Covered area
            area = cv2.contourArea(contours[0])*100/(img.shape[0]*img.shape[1])
            contours = np.vstack(contours).squeeze()
            self.statusText.set("Coverage: "+str(np.round(area*self.tv_res.get()**2*1e-6,8))+" mm\u00b2 ("+str(area)+" %)")
            
            mask = np.zeros(img.shape,np.uint8)
            
            #draw outline
            cv2.drawContours(mask,[contours],0,255,1)
            pixels = np.transpose(np.nonzero(mask))
            for i,p in enumerate(pixels):
                self.canvas_heat.create_oval(p[1], p[0], p[1], p[0], width = 1, fill="#007000", tags="o"+str(i))
                self.contourlist.append("o"+str(i))

        except Exception as e:
            print("error:", e)
        
    createVariables(self.master)
    createCanvasBinding()
    
def clearMarkingData(self):
    """
        Clear marking rectangle, outlines, and generated text
    """
    try:
        self.canvas_heat.delete('rectangleGA')
    except: 
        pass
    
    try:
        for i in self.contourlist:
            self.canvas_heat.delete(i)
        self.contourlist.clear()
    except: 
        pass