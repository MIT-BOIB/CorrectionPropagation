""" 
FrameManualCorrection: Module to create manual correction frame 
-------------------------------------------------------------------
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
from skimage import io
from tkinter import *
import PIL.Image, PIL.ImageTk
from skimage.transform import rescale as rsc
from Algorithms.Flattening.VolumeFlattening import unFlatten
from Algorithms.ManualCorrection.GraphCutBM import propagateBM
from Algorithms.ManualCorrection.GraphCutRPE import propagateRPE
from Algorithms.ManualCorrection.GraphCutILM import propagateILM
from Algorithms.ManualCorrection.Evaluation import calc_MSE_STDDEV
from matplotlib.colors import LinearSegmentedColormap
from FileHandler import ImportHandler
import threading
from Tooltips.TipHandler import CreateToolTip
import os
#threshold for mode
MODE_THRESHOLD = 16

def createManualCorrectionFrame(self):
    """ 
        Creation of Frame for manual refinement.
        
        This is the basic module for correction integration.
        If a new layer is added, add the abbreviation to the option menu at 
        the bottom (search: IntegrateOption).
        Furthermore, Algorithms -> ManualCorrection add a module 
        GraphCutxxx for the specific layer and add a call to propagate this 
        where 'IntegratePropagation' is commented. 
        You will also have to add a Layer segmentation value at line 33 and
        add segmentations where the comment 'IntegrateLayerVariable' are 
        located.
    
    Parameters:
    ----------
    self: object
        Framework
    """
    #diverging colormap
    c = ["blue","royalblue","white","indianred","red"]
    #fixed variables
    self.point_list = {}
    self.startpoint_x = None
    self.endpoint_x = None
    self.startpoint_y = None
    self.endpoint_y = None
    self.rectid = None
    self.segmentation_points = None
    self.savedslices = {}
    self._lines = []
    
    self._job_xz_corr = None
    self._job_volume_xy_corr = None
    self._job_HM = None
    
    self.parent = None
    self.rectx0 = 0
    self.recty0 = 0
    self.rectx1 = 0
    self.recty1 = 0
    self.move = False
    self.GT_EXISTS = True
    self.enface_lines = []
    #default segmentation layer
    self.layerSelection = "BM"
    
    # try to load ground truth automatically
    try:
        self.groundtruth = io.imread(self.initialdir+"//groundtruth.tif").astype('uint8')
        self.groundtruth = unFlatten(self.groundtruth, -self.shiftedValues).astype('uint8')
    except:
        self.GT_EXISTS = False
    
    def background(func):
        """ 
            Background worker - doesn't stop the GUI
            
            Parameters
            ----------
            func: function
                function to put in extra thread
        """
        th = threading.Thread(target=func)
        th.start()
    
    def layerChange(selection):
        """
            Change and inpaint correction layer.
            
            Parameters
            ----------
            selection: event, StringVar
                selected layer from dropdown menu
        """
        self.statusText.set("Change to:"+selection)
        self.layerSelection = selection
        
        resetVars()
        updateVolumeXZCorrSlice(None)
        
        if self.GT_EXISTS is True:
            createHeatmap()
        self.statusText.set("Change successful!")
    
    def saveSegmentation():
        """ 
            Save segmentation to directory. Saves as tif if nothing specified.
        """
        self.statusText.set("Saving Segmentation...")
        f =filedialog.asksaveasfilename(initialdir = self.initialdir,title = "Save as ...",filetypes = (("all files","*.*"),("tiff","*.tiff"),("tif","*.tif")))
        if f is None: 
            return
    
        if np.logical_not(f[-4:] == '.tif' or f[-5:] =='.tiff'):
            f += '.tif'

        saveVolume = self.segmentation.copy()
        try:
            saveVolume = unFlatten(saveVolume, self.shiftedValues).astype('uint8')
            self.statusText.set("Segmentation saved!")
        except Exception as e:
            self.statusText.set("Attention! Segmentation not saved!")
            print('Error:',e)
        io.imsave(f, saveVolume)
    
    def leftclickXYCorr(event):
        """ 
            Left click handler - Inpaint corrected line and store values
            
            Parameters
            ----------
            event: event
                click event
        """
        canvas_event = event.widget
        c_x = int(canvas_event.canvasx(event.x))
        c_y = int(canvas_event.canvasy(event.y))
        
        #dictionary entry
        self.savedslices[self.vol_sliderXZCorrection.get()] = np.zeros((1))
        
        if self.SMALLWINDOW is True:
            mul_fac = 2
            if(c_x % 2 == 1):
                c_x = (c_x - 1)//mul_fac
            else:
                c_x = c_x//mul_fac
        else:
            mul_fac = 3
            if(c_x % 3 == 1):
                c_x = (c_x - 1)//mul_fac
            elif(c_x % 3 == 2):
                c_x = (c_x - 2)//mul_fac
            else:
                c_x = c_x//mul_fac

        #inpaint new point and delete old - windows is strecthed by 3
        self.canvasXZCorrection.delete(str(c_x)+'a')
        self.canvasXZCorrection.delete(str(c_x)+'b')
        
        if mul_fac == 3:
            self.canvasXZCorrection.delete(str(c_x)+'c')

        self.segmentation_points[1][c_x] = c_y//mul_fac
        self.point_list[c_x]   = self.canvasXZCorrection.create_oval(c_x*mul_fac,   c_y, c_x*mul_fac,   c_y, width = 0, fill=self.crosshair_color, tags=str(c_x)+'a')
        self.point_list[c_x+1] = self.canvasXZCorrection.create_oval(c_x*mul_fac+1, c_y, c_x*mul_fac+1, c_y, width = 0, fill=self.crosshair_color, tags=str(c_x)+'b')
        if mul_fac == 3:
            self.point_list[c_x+2] = self.canvasXZCorrection.create_oval(c_x*mul_fac+2, c_y, c_x*mul_fac+2, c_y, width = 0, fill=self.crosshair_color, tags=str(c_x)+'c')
        
        #write into segmentation
        self.segmentation[self.vol_sliderXZCorrection.get(),:,c_x] = np.zeros((self.segmentation.shape[1]))
        # get mode, if wrong mode, BM
        ''' 
            IntegrateLayerVariable: Insert line when integrating new layer
        '''
        if self.layerSelection is 'BM':
            self.segmentation[self.vol_sliderXZCorrection.get(),c_y//mul_fac,c_x] = self.BM_VALUE
        elif self.layerSelection is 'RPE':
            self.segmentation[self.vol_sliderXZCorrection.get(),c_y//mul_fac,c_x] = self.RPE_VALUE
        elif self.layerSelection is 'ILM':
            self.segmentation[self.vol_sliderXZCorrection.get(),c_y//mul_fac,c_x] = self.ILM_VALUE
        else:
            self.segmentation[self.vol_sliderXZCorrection.get(),c_y//mul_fac,c_x] = self.BM_VALUE
            
        if self.startpoint_x is None or self.endpoint_x is None:
            self.startpoint_x = self.endpoint_x = c_x
        else:
            if c_x < self.startpoint_x:
                self.startpoint_x = c_x
            elif c_x > self.endpoint_x:
                self.endpoint_x = c_x
            if self.SMALLWINDOW is True:
                self.canvasXYCorrection.create_line(self.startpoint_x//2,self.vol_sliderXZCorrection.get()//2,self.endpoint_x//2,self.vol_sliderXZCorrection.get()//2,fill=self.crosshair_color,tags='enfaceline'+str(self.vol_sliderXZCorrection.get()))
            else:
                self.canvasXYCorrection.create_line(self.startpoint_x,self.vol_sliderXZCorrection.get(),self.endpoint_x,self.vol_sliderXZCorrection.get(),fill=self.crosshair_color,tags='enfaceline'+str(self.vol_sliderXZCorrection.get()))
            self.enface_lines.append('enfaceline'+str(self.vol_sliderXZCorrection.get()))
     
    def leftclickXYCorrMovement(event):
        """ 
            Left click and move handler - Inpaint corrected line and store values
            
            Compare function leftclickXYCorr(event)
            
            Parameters
            ----------
            event: event
                click event
        """
        canvas_event = event.widget
        c_x = int(canvas_event.canvasx(event.x))
        c_y = int(canvas_event.canvasy(event.y))
        
        if self.SMALLWINDOW is True:
            mul_fac = 2
            if(c_x % 2 == 1):
                c_x = (c_x - 1)//mul_fac
            else:
                c_x = c_x//mul_fac           
        else:
            mul_fac = 3
        
            if(c_x % 3 == 1):
                c_x = (c_x - 1)//mul_fac
            elif(c_x % 3 == 2):
                c_x = (c_x - 2)//mul_fac
            else:
                c_x = c_x//mul_fac
            
        self.canvasXZCorrection.delete(str(c_x)+'a')
        self.canvasXZCorrection.delete(str(c_x)+'b')
        if mul_fac == 3:
            self.canvasXZCorrection.delete(str(c_x)+'c')
        
        self.savedslices[self.vol_sliderXZCorrection.get()] = np.zeros((1))
        
        self.segmentation_points[1][c_x] = c_y//mul_fac
        self.point_list[c_x]   = self.canvasXZCorrection.create_oval(c_x*mul_fac,   c_y, c_x*mul_fac,   c_y, width = 0, fill=self.crosshair_color, tags=str(c_x)+'a')
        self.point_list[c_x+1] = self.canvasXZCorrection.create_oval(c_x*mul_fac+1, c_y, c_x*mul_fac+1, c_y, width = 0, fill=self.crosshair_color, tags=str(c_x)+'b')
        if mul_fac == 3:
            self.point_list[c_x+2] = self.canvasXZCorrection.create_oval(c_x*mul_fac+2, c_y, c_x*mul_fac+2, c_y, width = 0, fill=self.crosshair_color, tags=str(c_x)+'c')
                
        self.segmentation[self.vol_sliderXZCorrection.get(),:,c_x] = np.zeros((self.segmentation.shape[1]))
        ''' 
            IntegrateLayerVariable: Insert line when integrating new layer
        '''
        if self.layerSelection is 'BM':
            self.segmentation[self.vol_sliderXZCorrection.get(),c_y//mul_fac,c_x] = self.BM_VALUE
        elif self.layerSelection is 'RPE':
            self.segmentation[self.vol_sliderXZCorrection.get(),c_y//mul_fac,c_x] = self.RPE_VALUE
        elif self.layerSelection is 'ILM':
            self.segmentation[self.vol_sliderXZCorrection.get(),c_y//mul_fac,c_x] = self.ILM_VALUE
        else:
            self.segmentation[self.vol_sliderXZCorrection.get(),c_y//mul_fac,c_x] = self.BM_VALUE
        
        if self.startpoint_x is None or self.endpoint_x is None:
            self.startpoint_x = self.endpoint_x = c_x
        else:
            if c_x < self.startpoint_x:
                self.startpoint_x = c_x
            elif c_x > self.endpoint_x:
                self.endpoint_x = c_x
            if self.SMALLWINDOW is True:
                self.canvasXYCorrection.create_line(self.startpoint_x//2,self.vol_sliderXZCorrection.get()//2,self.endpoint_x//2,self.vol_sliderXZCorrection.get()//2,fill=self.crosshair_color,tags='enfaceline'+str(self.vol_sliderXZCorrection.get()))
            else:
                self.canvasXYCorrection.create_line(self.startpoint_x,self.vol_sliderXZCorrection.get(),self.endpoint_x,self.vol_sliderXZCorrection.get(),fill=self.crosshair_color,tags='enfaceline'+str(self.vol_sliderXZCorrection.get()))
            self.enface_lines.append('enfaceline'+str(self.vol_sliderXZCorrection.get()))
            
    def resetSlice():
        """ 
            Reset correction of current slice
        """
        self.statusText.set("Resetting slice...")
        try:
            self.segmentation[self.vol_sliderXZCorrection.get(),:,:] = self.segmentation_original[self.vol_sliderXZCorrection.get(),:,:]
            del self.savedslices[self.vol_sliderXZCorrection.get()]
            updateVolumeXZCorrSlice(None)
            
            try:
                self.canvasXYCorrection.delete('enfaceline'+str(self.vol_sliderXZCorrection.get()))
                self.enface_lines.remove('enfaceline'+str(self.vol_sliderXZCorrection.get()))
            except:
                pass
        except:
            pass
        self.statusText.set("Slice reset done.")
    
    def inpaintCorrSegmentation():
        """ 
            Inpaint corrected BM segmentation to volume
        """
        #select layer
        ''' 
            IntegrateLayerVariable: Insert line when integrating new layer
        '''
        if self.layerSelection is 'BM':
            coordinates = np.asarray(np.where(self.segmentation[self.vol_sliderXZCorrection.get()] == self.BM_VALUE))
        elif self.layerSelection is 'RPE':
            coordinates = np.asarray(np.where(self.segmentation[self.vol_sliderXZCorrection.get()] == self.RPE_VALUE))
        elif self.layerSelection is 'ILM':
            coordinates = np.asarray(np.where(self.segmentation[self.vol_sliderXZCorrection.get()] == self.ILM_VALUE))
        else:
            coordinates = np.asarray(np.where(self.segmentation[self.vol_sliderXZCorrection.get()] == self.BM_VALUE))
        
        if self.SMALLWINDOW is True:
            mul_fac = 2
        else:
            mul_fac = 3

        dict = {}
        for x in range(coordinates.shape[1]):
            dict[coordinates[1][x]] = coordinates[0][x]
        
        for key,value in dict.items():
            self.point_list[str(key)+'a'] = self.canvasXZCorrection.create_oval(key*mul_fac, value*mul_fac, key*mul_fac, value*mul_fac, width = 0, fill=self.crosshair_color, tags=str(key)+'a')
            self.point_list[str(key)+'b'] = self.canvasXZCorrection.create_oval(key*mul_fac+1, value*mul_fac, key*mul_fac+1, value*mul_fac, width = 0, fill=self.crosshair_color, tags=str(key)+'b')
            if mul_fac == 3:
                self.point_list[str(key)+'c'] = self.canvasXZCorrection.create_oval(key*mul_fac+2, value*mul_fac, key*mul_fac+2, value*mul_fac, width = 0, fill=self.crosshair_color, tags=str(key)+'c')
            
    def _createVariables(parent_var):
        """ 
            Create basic variables
        """
        self.parent = parent_var
        self.rectx0 = 0
        self.recty0 = 0
        self.rectx1 = 0
        self.recty1 = 0
        self.rectid = None
        self.move = False
    
    def _createCanvasBinding():
        """ 
            Create basic bindings
        """
        self.canvasXYCorrection.bind( "<Button-1>", startRect )
        self.canvasXYCorrection.bind( "<ButtonRelease-1>", stopRect )
        self.canvasXYCorrection.bind( "<Motion>", movingRect )

    def startRect(event):
        """ 
            start rectangular selection
            
            Parameters
            ----------
            event: event
                click event
        """
        self.move = True
        self.canvasXYCorrection.delete('enFaceRect') 
        self.rectx0 = self.canvasXYCorrection.canvasx(event.x)
        self.recty0 = self.canvasXYCorrection.canvasy(event.y) 

        self.rect = self.canvasXYCorrection.create_rectangle(self.rectx0, self.recty0, self.rectx0, self.recty0, outline=self.crosshair_color,width=3, tags='enFaceRect')
        self.rectid = self.canvasXYCorrection.find_closest(self.rectx0, self.recty0, halo=2)

    def movingRect(event):
        """ 
            Move endpoint of rectangle
            
            Parameters
            ----------
            event: event
                click event
        """
        if self.move == True: 
            self.rectx1 = self.canvasXYCorrection.canvasx(event.x)
            self.recty1 = self.canvasXYCorrection.canvasy(event.y)
            self.canvasXYCorrection.coords(self.rectid, self.rectx0, self.recty0, self.rectx1, self.recty1)

    def stopRect(event):
        """ 
            Stop rectangular selection and check if valid
            
           Parameters
            ----------
            event: event
                click event
        """
        try:
            self.move = False
            self.rectx1 = self.canvasXYCorrection.canvasx(event.x)
            self.recty1 = self.canvasXYCorrection.canvasy(event.y) 
            
            if self.SMALLWINDOW is True:
                self.canvasXYCorrection.coords(self.rectid, self.startpoint_x//2, self.recty0, self.endpoint_x//2, self.recty1)
                self.box = (int(self.startpoint_x), int(self.recty0), int(self.endpoint_x)-int(self.startpoint_x), int(self.recty1)-int(self.recty0))
            else:
                self.canvasXYCorrection.coords(self.rectid, self.startpoint_x, self.recty0, self.endpoint_x, self.recty1)
                self.box = (int(self.startpoint_x), int(self.recty0), int(self.endpoint_x)-int(self.startpoint_x), int(self.recty1)-int(self.recty0))
          
            if(checkRectValid() == False):
                self.canvasXYCorrection.delete('enFaceRect') 
                self.startpoint_y = self.endpoint_y  = self.vol_sliderXZCorrection.get()
            else:
                if self.SMALLWINDOW is True:
                    self.startpoint_y = int(self.recty0*2)  
                    self.endpoint_y = int(self.recty1*2)  
                else:
                    self.startpoint_y = int(self.recty0)  
                    self.endpoint_y = int(self.recty1) 
        except:
            pass
                        
    def checkRectValid():
        """ 
            Check if rectangle is valid. It has to cover at least one point on each line.
        """
        if self.SMALLWINDOW is True:
            if (self.recty1*2 >= self.vol_sliderXZCorrection.get() and self.recty0*2 <= self.vol_sliderXZCorrection.get()) or (self.recty1*2 <= self.vol_sliderXZCorrection.get() and self.recty0*2 >= self.vol_sliderXZCorrection.get()):
                return True
            else:
                return False
        else:
            if (self.recty1 >= self.vol_sliderXZCorrection.get() and self.recty0 <= self.vol_sliderXZCorrection.get()) or (self.recty1 <= self.vol_sliderXZCorrection.get() and self.recty0 >= self.vol_sliderXZCorrection.get()):
                return True
            else:
                return False
            
    def propagateCorrection():
        """ 
            Propagate manual corrections.
            
            The segmentation is cropped and based on the corrections, l and r are calculated. From this,
            the mode is determined and the linear approximation (l<16) or graph-cut (l>=16) executed.
        """
        try:
            
            #get corrected slices
            for key in self.savedslices:
                self.savedslices[key] = self.segmentation[key]
            #determine mode
            mode = 'high'
            if (self.endpoint_y-self.startpoint_y)//len(self.savedslices) < MODE_THRESHOLD:
                mode = 'low'
                
            self.statusText.set("Running: (r="+str(self.endpoint_x-self.startpoint_x)+'x'+str(self.endpoint_y-self.startpoint_y)+ ", l="+str((self.endpoint_y-self.startpoint_y)//len(self.savedslices))+")")
            #propagation algorithm
            '''
            IntegratePropagation: Add lines for new layers
            '''
            rect_correction = (self.startpoint_x,self.endpoint_x,self.startpoint_y,self.endpoint_y)
            
            if self.layerSelection is 'BM':
                cropped_segmentation = propagateBM(self.volume_original, self.segmentation, self.segmentation_original, self.savedslices, rect_correction, mode)
            elif self.layerSelection is 'RPE':
                cropped_segmentation = propagateRPE(self.volume_original, self.segmentation, self.segmentation_original, self.savedslices, rect_correction, mode)
            elif self.layerSelection is 'ILM':
                cropped_segmentation = propagateILM(self.volume_original, self.segmentation, self.segmentation_original, self.savedslices, rect_correction, mode)
            else:
                cropped_segmentation = propagateBM(self.volume_original, self.segmentation, self.segmentation_original, self.savedslices, rect_correction, mode)    
  
            self.segmentation[self.startpoint_y-1:self.endpoint_y+1,:,self.startpoint_x-1:self.endpoint_x+2] = cropped_segmentation
            #update slices
            updateVolumeXZCorrSlice(None)
            updateCorrectionHeatmap()
            
            resetVars()

        except Exception as e:
            print("Error:",e)
            updateVolumeXZCorrSlice(None)
            updateCorrectionHeatmap()
            resetVars()
            
        self.statusText.set("Propagation finished!")
        
            
    def createHeatmap():
        """ 
            Create canvas when ground truth loaded/exists
            
            The colormap is diverging, red/blue are poistive/negative distances, whereas white is a perfect match.
            This tool can be used to evaluate an automatic segmentation algorithm.
        """
        if hasattr(self, "FrameEvaluation"):
            self.FrameEvaluation.destroy()
            del self.FrameEvaluation
            
        self._job_HM = None
        
        if self.SMALLWINDOW is True:
            self.FrameEvaluation = Frame(self.correctionFrame,width = self.widthXY//2+30, height = self.heightXY//2,bg=self.canvasbackground)
            self.canvasEvaluation = Canvas(self.FrameEvaluation, width = self.widthXY//2, height = self.heightXY//2,highlightthickness=2,highlightbackground=self.hlghtbg,bg=self.canvasbackground)    
            self.canvasColorbarCorr = Canvas(self.FrameEvaluation, width = 1, height = self.heightXY//2,bg=self.canvasbackground)
            self.canvas_Corrcbar_HM = Canvas(self.canvasColorbarCorr, width=20, height = self.heightXY//2 - 40,bg=self.canvasbackground)
            self.btnAutoEval = Button(self.FrameEvaluation, width = 18 , height = 2, text="Automatic Evaluation",font=self.statusFont, command=lambda:background(propagateCorrectionEvaluation), activebackground=self.btn_common_bg, fg=self.btnforeground, bg=self.btnbackground, activeforeground=self.btnforeground)
            self.btnAutoEval.grid(row=0,column=0,padx=(10,10),pady=(215,0), sticky="NW")
        else:
            self.FrameEvaluation = Frame(self.correctionFrame,width = self.widthXY+30, height = self.heightXY,bg=self.canvasbackground)
            self.canvasEvaluation = Canvas(self.FrameEvaluation, width = self.widthXY, height = self.heightXY,highlightthickness=2,highlightbackground=self.hlghtbg,bg=self.canvasbackground)
            self.canvasColorbarCorr = Canvas(self.FrameEvaluation, width = 1, height = self.heightXY,bg=self.canvasbackground)
            self.canvas_Corrcbar_HM = Canvas(self.canvasColorbarCorr, width=20, height = self.heightXY - 40,bg=self.canvasbackground)
            self.btnAutoEval = Button(self.FrameEvaluation, width = 15 , height = 3, text="Automatic Evaluation",font=self.statusFont, command=lambda:background(propagateCorrectionEvaluation), activebackground=self.btn_common_bg, fg=self.btnforeground, bg=self.btnbackground, activeforeground=self.btnforeground)
            self.btnAutoEval.grid(row=0,column=0,padx=(10,10),pady=(435,0), sticky="NW")
        
        self.FrameEvaluation.grid(row=0,column=2, sticky="NE")
        
        self.canvasEvaluation.grid(row=0,column=1, sticky="W")
        CreateToolTip(self.btnAutoEval, self.ttip_dict['manAutoEval'])
        
        #create difference map
        self.differencemap = np.zeros((self.groundtruth.shape[0],self.groundtruth.shape[2])).astype('int32')
        ''' 
            IntegrateLayerVariable: Insert line when integrating new layer
        '''        
        for z in range(self.groundtruth.shape[0]):
            for x in range (self.groundtruth.shape[2]):
                if self.layerSelection is 'BM':
                    values_algorithm  = np.where(self.segmentation[z,:,x] == self.BM_VALUE)[0]
                    values_gt         = np.where(self.groundtruth[z,:,x] == self.BM_VALUE)[0]
                elif self.layerSelection is 'RPE':
                    values_algorithm  = np.where(self.segmentation[z,:,x] == self.RPE_VALUE)[0]
                    values_gt         = np.where(self.groundtruth[z,:,x] == self.RPE_VALUE)[0]
                elif self.layerSelection is 'ILM':
                    values_algorithm  = np.where(self.segmentation[z,:,x] == self.ILM_VALUE)[0]
                    values_gt         = np.where(self.groundtruth[z,:,x] == self.ILM_VALUE)[0]
                else:
                    values_algorithm  = np.where(self.segmentation[z,:,x] == self.BM_VALUE)[0]
                    values_gt         = np.where(self.groundtruth[z,:,x] == self.BM_VALUE)[0]

                try:
                    self.differencemap[z,x] = values_gt - values_algorithm
                except:
                    self.differencemap[z,x] = 0

        if(np.max(self.differencemap) <= 0):
            range_of_values = np.abs(np.max(self.differencemap))+np.abs(np.min(self.differencemap))
        else:
            range_of_values = np.abs(np.max(self.differencemap))+np.abs(np.min(self.differencemap)) + 1
        
        offset = np.abs(np.min(self.differencemap))
        
        if (range_of_values == 0): 
            range_of_values = 8
            offset = 4
            
        self.differencemap += offset
        self.differencemap = self.differencemap.astype('uint8')
        divisor = float(offset/range_of_values)
        
        #colorbar
        v = [0,divisor/2,divisor, divisor + (1-divisor)/2,1.]
        l = list(zip(v,c))
        
        if self.SMALLWINDOW is True:
            self.differencemap = self.differencemap[::2,::2]
        
        self.cmap_corr = LinearSegmentedColormap.from_list('rg',l, N=range_of_values)
        self.heatmap_diff = PIL.ImageTk.PhotoImage(image = PIL.Image.fromarray(self.cmap_corr(self.differencemap, bytes=True)))
        self.image_on_canvas_diff = self.canvasEvaluation.create_image(0, 0, image=self.heatmap_diff, anchor=NW)
        
        self.canvasColorbarCorr.grid(row=0,column=2,sticky="nw")
        
        self.label_Corrtop_HM = Label(self.canvasColorbarCorr, width = 3,height= 1, text=''+str(int(range_of_values-offset)),font=self.clrbarFont,anchor="nw",bg=self.canvasbackground, fg =self.canvasforeground)
        self.label_Corrtop_HM.grid(row=0,column=0)
        self.canvas_Corrcbar_HM.grid(row=1,column=0,padx=(3, 0),pady = (0,0),sticky="nw")
        self.label_Corrbottom_HM = Label(self.canvasColorbarCorr, width = 3, height=  1, text=""+str(int(-offset)),font=self.clrbarFont,anchor="nw",bg=self.canvasbackground, fg =self.canvasforeground)
        self.label_Corrbottom_HM.grid(row=2,column=0)
        
        cbar_arr = np.zeros((self.heightXY, 20)).astype('uint8')
        incrementer = range_of_values/self.heightXY
        for i in range (self.heightXY):
            cbar_arr[i,:] = incrementer*i 
        #update canvas
        if self.SMALLWINDOW is True:
            cbar_arr = cbar_arr[::2,:]
            
        self.photo_cbar_corr = PIL.ImageTk.PhotoImage(image = PIL.Image.fromarray(self.cmap_corr(cbar_arr[::-1], bytes=True)))
        self.image_on_canvasCorr_cbar = self.canvas_Corrcbar_HM.create_image(0, 0, image=self.photo_cbar_corr, anchor=NW)
        
    def updateVolumeXZCorrSlice(event=None): 
        """ 
            Update XZ correction slice (time delay)
            
            Optional
            ---------
            event: event
                event (unused)
        """
        if self._job_xz_corr:
            self.master.after_cancel(self._job_xz_corr)
        self._job_xz_corr = self.master.after(50, updateVolumeXZCorrSlice_Helper)
        
    def updateVolumeXZCorrSlice_Helper():
        """ 
            Update XZ correction slice helper
            
            Load new slice, resize by factor 3 and paint to canvas.
        """
        self.canvasXZCorrection.delete("all")

        self._job_xz_corr = None
        
        self.segmentation_slice  = self.segmentation[self.vol_sliderXZCorrection.get()].astype('uint8')
        ''' 
            IntegrateLayerVariable: Insert line when integrating new layer
        '''
        if self.layerSelection is 'BM':
            self.segmentation_points = np.asarray(np.where(self.segmentation_slice.transpose() == self.BM_VALUE))
        elif self.layerSelection is 'RPE':
            self.segmentation_points = np.asarray(np.where(self.segmentation_slice.transpose() == self.RPE_VALUE))
        elif self.layerSelection is 'ILM':
            self.segmentation_points = np.asarray(np.where(self.segmentation_slice.transpose() == self.ILM_VALUE))
        else:
            self.segmentation_points = np.asarray(np.where(self.segmentation_slice.transpose() == self.BM_VALUE))
        
        self.sliceXZCorrection = self.volumeXZ[self.vol_sliderXZCorrection.get()].astype('uint8')
        if self.SMALLWINDOW is True:
            self.sliceXZCorrection = rsc(self.sliceXZCorrection , 2, preserve_range = True).astype('uint8')
        else:
            self.sliceXZCorrection = rsc(self.sliceXZCorrection , 3, preserve_range = True).astype('uint8')
        self.slice_pilXZCorrection = PIL.Image.fromarray(self.sliceXZCorrection)  
        
        self.photoXZCorrection = PIL.ImageTk.PhotoImage(image = self.slice_pilXZCorrection)
        self.image_on_canvas_XZCorrection = self.canvasXZCorrection.create_image(0, 0, image=self.photoXZCorrection, anchor=NW)
        self.canvasXZCorrection.configure(scrollregion = self.canvasXZCorrection.bbox("all"))
        inpaintCorrSegmentation()
        
    def updateVolumeSliceCorr(event):
        """ 
            Update en face correction slice (time delay)
            
            Optional
            ---------
            event: event
                event (unused)
        """
        if self._job_volume_xy_corr :
            self.master.after_cancel(self._job_volume_xy_corr)
        self._job_volume_xy_corr  = self.master.after(50, updateVolumeSliceCorr_Helper)
            
    def updateVolumeSliceCorr_Helper():
        """ 
            Update en face correction slice helper
            
            Load new slice
        """
        self._job_volume_xy_corr  = None

        self.sliceEnface = self.volume[self.vol_sliderXY_corr.get()].astype('uint8')
        
        if self.SMALLWINDOW is True:
            self.sliceEnface = self.sliceEnface[::2,::2]
            
        self.slice_pil_Enface = PIL.Image.fromarray(self.sliceEnface)   
        self.photoEnface = PIL.ImageTk.PhotoImage(image = self.slice_pil_Enface)
        self.canvasXYCorrection.itemconfig(self.image_on_canvas_Enface, image = self.photoEnface)
        self.canvasXYCorrection.configure(scrollregion = self.canvasXYCorrection.bbox("all"))
        
        if self.mover_init is True:
            if self.SMALLWINDOW is True:
                self.canvasXZCorrection.yview_moveto((self.heightXZ*2//2-175)/(self.heightXZ*2))
            else:
                self.canvasXZCorrection.yview_moveto((self.heightXZ*3//2-175)/(self.heightXZ*3))
            self.mover_init = False
            
    def updateCorrectionHeatmap():
        """ 
            Update ground truth heatmap slice (time delay)
        """
        if self._job_HM:
            self.master.after_cancel(self._job_HM)
        self._job_HM = self.master.after(50, updateCorrectionHeatmap_Helper)
        
    def updateCorrectionHeatmap_Helper():
        """ 
            Update ground truth heatmap slice helper
            
            Load new slice
        """
        self._job_HM = None

        self.differencemap = np.zeros((self.groundtruth.shape[0],self.groundtruth.shape[2])).astype('int32')
        ''' 
            IntegrateLayerVariable: Insert line when integrating new layer
        '''     
        for z in range(self.groundtruth.shape[0]):
            for x in range (self.groundtruth.shape[2]):
                if self.layerSelection is 'BM':
                    values_algorithm  = np.where(self.segmentation[z,:,x] == self.BM_VALUE)[0]
                    values_gt         = np.where(self.groundtruth[z,:,x] == self.BM_VALUE)[0]
                elif self.layerSelection is 'RPE':
                    values_algorithm  = np.where(self.segmentation[z,:,x] == self.RPE_VALUE)[0]
                    values_gt         = np.where(self.groundtruth[z,:,x] == self.RPE_VALUE)[0]
                elif self.layerSelection is 'ILM':
                    values_algorithm  = np.where(self.segmentation[z,:,x] == self.ILM_VALUE)[0]
                    values_gt         = np.where(self.groundtruth[z,:,x] == self.ILM_VALUE)[0]
                else:
                    values_algorithm  = np.where(self.segmentation[z,:,x] == self.BM_VALUE)[0]
                    values_gt         = np.where(self.groundtruth[z,:,x] == self.BM_VALUE)[0]

                try:
                    self.differencemap[z,x] = values_gt - values_algorithm
                except:
                    self.differencemap[z,x] = 0
        
        offset = np.abs(np.min(self.differencemap))
        self.differencemap += offset
        self.differencemap = self.differencemap.astype('uint8')
        
        if self.SMALLWINDOW is True:
            self.differencemap = self.differencemap[::2,::2]

        self.heatmap_diff = PIL.ImageTk.PhotoImage(image = PIL.Image.fromarray(self.cmap_corr(self.differencemap, bytes=True)))
        
        self.canvasEvaluation.itemconfig(self.image_on_canvas_diff, image = self.heatmap_diff)
    
    def loadGroundtruth():
        """ 
            Load ground truth volume
        """
        self.statusText.set("Loading GT...")
        try:
            self.groundtruth = ImportHandler.loadGroundtruthVolume(self.initialdir,self.shiftedValues)
            self.GT_EXISTS = True
            if(hasattr(self, "FrameEvaluation")):
                del self.FrameEvaluation
            createHeatmap()
            self.statusText.set("GT loaded!")
        except Exception as e :
            self.statusText.set("Error when loading GT!")
            print("Error:",e)
            self.GT_EXISTS = False
        
    def propagateCorrectionEvaluation():
        """ 
            Run an automatic algorithm evaluation of the current volume.
            
            Define r,l, and spacing in auto_evaluation.txt. The evaluation then 
            evaluates MSE and STDDEV for the algorithm of the current volume. The 
            algorithm uses the ground truth as input for the corrected lines and
            sets an equal spacing.
            
            Example:
            1) Load Volume
            2) Run automatic segmentation
            3) Switch to Manual Refinement
            4) Load ground truth
            
            Parameters in txt:
            r=32,128
            l=2,16
            spacing=32
            
            5)  An evaluation of the manual correction algorithm is run for the 
                following settings:
                
                Two different region sizes (r): 32x32 and 128x128
                In those regions, two different line settings are evaluated (l): 
                    i) every second line 
                    ii) every 16th line
                is corrected by inserting the underlying ground truth lines
            
            6) Spacing describes the stride from of the batches. The lower the number,
               the higher the computation time. For initial experiments, half of the 
               region r is sufficient!!!             
            
        """
        self.statusText.set("Running evaluation...")
        s = ttk.Style()
        s.theme_use('clam')
        s.configure("PRLEC.Horizontal.TProgressbar", troughcolor='#A0A0A0', background="#32E3E3", thickness=1, troughrelief="flat", relief="flat", borderwidth=0)
        
        try:
            if(np.count_nonzero(self.groundtruth) == 0):
                return
            f = open(os.path.dirname(os.path.abspath(sys.argv[0]))+"\\auto_evaluation.txt", "r")
            dimensions = f.readline().replace("r=","").replace("\n","").split(',')
            correction_lines = f.readline().replace("l=","").replace("\n","").split(',')
            spacing = f.readline().replace("spacing=","").replace("\n","")
            f.close()
        except Exception as e:
            print("Error when running automatic evaluation:\n",e)
            return
        
        dim_val=IntVar()
        dim_val.set(1)
        self.labelprogressBarDimension = Label(self.FrameEvaluation, width = 12,height= 1, text='Region progress',font=self.textFont,anchor="nw",bg=self.canvasbackground, fg =self.canvasforeground) 
        self.progressBarDimension = ttk.Progressbar(self.FrameEvaluation, style="PRLEC.Horizontal.TProgressbar",orient="horizontal",length=140, maximum=len(dimensions), mode='determinate', variable = dim_val)
        
        if self.SMALLWINDOW is True:
            self.labelprogressBarDimension.grid(row=0,column=0,padx=(10,10),pady=(80,0), sticky="NW")
            self.progressBarDimension.grid(row=0,column=0,padx=(10,10),pady=(110,0), sticky="NW")
        else:
            self.labelprogressBarDimension.grid(row=0,column=0,padx=(10,10),pady=(280,0), sticky="NW")
            self.progressBarDimension.grid(row=0,column=0,padx=(10,10),pady=(320,0), sticky="NW")            
        
        f = open(os.path.dirname(os.path.abspath(sys.argv[0]))+"\\Evaluation.txt", "a")
        ''' 
            IntegrateLayerVariable: Insert line when integrating new layer
        '''
        if self.layerSelection is 'BM':
            seg_val = self.BM_VALUE
            f.write("Bruchs membrane\n")
        elif self.layerSelection is 'RPE':
            seg_val = self.RPE_VALUE
            f.write("Retinal pigment epithelium\n")
        elif self.layerSelection is 'ILM':
            seg_val = self.ILM_VALUE
            f.write("Inner limiting membrane\n")
        else:
            seg_val = self.BM_VALUE
            f.write("Bruchs membrane\n")
            
        mean_val,stddev_val = calc_MSE_STDDEV(self.segmentation, self.groundtruth, seg_val)
        f.write(self.initialdir+"\n")
        f.write("FA MSE: "+str(mean_val)+" Stddev: "+str(stddev_val)+"\n")
        spacing = int(spacing)
        
        self.labelprogressBarLines = Label(self.FrameEvaluation, width = 12,height= 1, text='Line progress',font=self.textFont,anchor="nw",bg=self.canvasbackground, fg =self.canvasforeground)

        region_val=IntVar()
        region_val.set(1)
        self.progressBarLines = ttk.Progressbar(self.FrameEvaluation, style="PRLEC.Horizontal.TProgressbar", orient="horizontal", length=140, maximum=100, mode='determinate', variable=region_val)

        
        if self.SMALLWINDOW is True:
            self.labelprogressBarLines.grid(row=0,column=0,padx=(10,10),pady=(150,0), sticky="NW")
            self.progressBarLines.grid(row=0,column=0,padx=(10,10),pady=(180,0), sticky="NW")
        else:
            self.labelprogressBarLines.grid(row=0,column=0,padx=(10,10),pady=(360,0), sticky="NW")
            self.progressBarLines.grid(row=0,column=0,padx=(10,10),pady=(400,0), sticky="NW")
        
        for dim_ctr,dim in enumerate(dimensions):
            dim_val.set(dim_ctr+1)
            dim = int(dim)            
            
            gt_mean_list = []
            gt_stddev_list = []
            
            test_mean_list = []
            test_stddev_list = []
            
            for no_of_l in range(len(correction_lines)):
                test_mean_list.append([])
                test_stddev_list.append([])
    
            f.write("r="+str(dim)+"\n")

            self.progressBarLines.configure(maximum = int((self.volume.shape[2] - dim)/spacing * (self.volume.shape[1] - dim - 1)/ spacing)+1)
            ctr_lines = 1
            for x in range(1, self.volume.shape[2] - dim ,spacing):
                
                self.startpoint_x =  x
                self.endpoint_x = self.startpoint_x + dim - 1
                
                for y in range(1, self.volume.shape[1] - dim - 1, spacing):
                    region_val.set(ctr_lines)
                    self.startpoint_y =  y
                    self.endpoint_y = self.startpoint_y + dim
         
                    cropped_segmentation = self.segmentation[self.startpoint_y-1:self.endpoint_y+1,:,self.startpoint_x-1:self.endpoint_x+2].astype('uint8')
                    cropped_segmentation[0,:,:] = self.groundtruth[self.startpoint_y-1,:,self.startpoint_x-1:self.endpoint_x+2]
                    cropped_segmentation[-1,:,:] = self.groundtruth[self.endpoint_y+1,:,self.startpoint_x-1:self.endpoint_x+2]
                    
                    ground_truth_crop = self.groundtruth[self.startpoint_y:self.endpoint_y,:,self.startpoint_x:self.endpoint_x+1]
    
                    mean_val,stddev_val = calc_MSE_STDDEV(cropped_segmentation[1:-1,:,1:-1], ground_truth_crop,seg_val)
    
                    gt_mean_list.append(mean_val)
                    gt_stddev_list.append(stddev_val)
    
                    for ctr,i in enumerate(correction_lines):
                        i = int(i)
                        mode ='low'
                        if(i > 16):
                            mode = 'high'
                        ''' 
                            IntegrateLayerVariable: Insert line when integrating new layer
                        '''                           
                        for z in range(self.startpoint_y+i,self.endpoint_y+1, i):
                            if self.layerSelection is 'BM':
                                self.savedslices[z] = np.where(self.groundtruth[z] == self.BM_VALUE, self.BM_VALUE,0)
                            elif self.layerSelection is 'RPE':
                                self.savedslices[z] = np.where(self.groundtruth[z] == self.RPE_VALUE, self.RPE_VALUE,0) 
                            elif self.layerSelection is 'ILM':
                                self.savedslices[z] = np.where(self.groundtruth[z] == self.ILM_VALUE, self.ILM_VALUE,0)
                            else:
                                self.savedslices[z] = np.where(self.groundtruth[z] == self.BM_VALUE, self.BM_VALUE,0)
                            
                        rect_correction = (self.startpoint_x,self.endpoint_x,self.startpoint_y,self.endpoint_y)
                        '''
                        IntegratePropagation: Add lines for new layers
                        ''' 
                        if self.layerSelection is 'BM':
                            cropped_segmentation = propagateBM(self.volume_original, self.segmentation, self.segmentation_original, self.savedslices, rect_correction, mode)
                        elif self.layerSelection is 'RPE':
                            cropped_segmentation = propagateRPE(self.volume_original, self.segmentation, self.segmentation_original, self.savedslices, rect_correction, mode)
                        elif self.layerSelection is 'ILM':
                            cropped_segmentation = propagateILM(self.volume_original, self.segmentation, self.segmentation_original, self.savedslices, rect_correction, mode)
                        else:
                            cropped_segmentation = propagateBM(self.volume_original, self.segmentation, self.segmentation_original, self.savedslices, rect_correction, mode)                        
                            
                        mean,stddev = calc_MSE_STDDEV(cropped_segmentation[1:-1,:,1:-1], ground_truth_crop,seg_val)
                        
                        test_mean_list[ctr].append(mean)
                        test_stddev_list[ctr].append(stddev)
                        
                        self.savedslices.clear()
                    ctr_lines += 1
            
            f.write("AveragePatch MSE: "+str(np.mean(np.array(gt_mean_list)))+" SDEV: " + str(np.mean(np.array(gt_stddev_list)))+'\n')
            for no_of_l in range(len(correction_lines)):
                f.write("L="+correction_lines[no_of_l]+" MSE: "+str(np.mean(np.array(test_mean_list[no_of_l]))) + " SDEV: "+str(np.mean(np.array(test_stddev_list[no_of_l])))+'\n')  
        
        self.progressBarLines.destroy()
        self.progressBarDimension.destroy()
        self.labelprogressBarLines.destroy()
        self.labelprogressBarDimension.destroy()
        
        f.close()
        self.statusText.set("")
        
    def resetVars():
        """
            Reset all variables necessary to correct
        """
        self.point_list = {}
        self.startpoint_x = None
        self.endpoint_x = None
        self.startpoint_y = None
        self.endpoint_y = None
        self.rectid = None
        self.segmentation_points = None
        self.savedslices = {}
        #delete correction lines
        for line in self.enface_lines:
            self.canvasXYCorrection.delete(line)
        self.enface_lines = []
            
        self.canvasXYCorrection.delete('enFaceRect')
    
    '''
        Basic frame setup
    '''
    if self.SMALLWINDOW is True:       
        self.correctionFrame = Frame(self.master, width = self.widthXZ, height = self.heightXY//2 + self.heightXZ + 25, relief=SUNKEN,bg=self.canvasbackground)
        self.canvasXYCorrection = Canvas(self.correctionFrame, width = self.widthXY//2+5, height = self.heightXY//2,highlightthickness=2,highlightbackground=self.hlghtbg,bg=self.canvasbackground)
        self.sliceEnface = self.volume[self.volume.shape[0]//2][::2,::2].astype('uint8')
        self.vol_sliderXY_corr = Scale(self.correctionFrame, from_=0, to=self.volume.shape[0]-1, length =self.widthXZ//2-4, orient=VERTICAL, command=updateVolumeSliceCorr,bg=self.canvasbackground, fg =self.canvasforeground)
        pad_x_sliderXY = (475,0)
        self.canvasXZCorrection = Canvas(self.correctionFrame, width = self.widthXZ*2,   height = 350,bg=self.canvasbackground,highlightthickness=2,highlightbackground=self.hlghtbg)
        self.manCorrBtnHeight = 1
        
        self.corr_y_pad_increment = 40
    else:
        self.correctionFrame = Frame(self.master, width = self.widthXZ*3, height = self.heightXY + self.heightXZ*3 + 50, relief=SUNKEN,bg=self.canvasbackground)
        self.canvasXYCorrection = Canvas(self.correctionFrame, width = self.widthXY+5, height = self.heightXY,highlightthickness=2,highlightbackground=self.hlghtbg,bg=self.canvasbackground)
        self.sliceEnface = self.volume[self.volume.shape[0]//2].astype('uint8')
        self.vol_sliderXY_corr = Scale(self.correctionFrame, from_=0, to=self.volume.shape[0]-1, length =self.widthXZ-4, orient=VERTICAL, command=updateVolumeSliceCorr,bg=self.canvasbackground, fg =self.canvasforeground)
        pad_x_sliderXY = (710,0)
        self.canvasXZCorrection = Canvas(self.correctionFrame, width = self.widthXZ*3, height = 350,bg=self.canvasbackground,highlightthickness=2,highlightbackground=self.hlghtbg)
        self.manCorrBtnHeight = 2
        self.corr_y_pad_increment = 80
    
    self.correctionFrame.grid(row = 0, column = 2, sticky ='nw', pady=(10,0))
    #setting bruchs correction as default
    self.varManCorrLayer = StringVar(self.correctionFrame)
    self.varManCorrLayer.set("BM")

    self.canvasXYCorrection.grid(row=0,column=0, padx=(200,0), sticky="NW")

    self.slice_pil_Enface = PIL.Image.fromarray(self.sliceEnface)
    self.photoEnface = PIL.ImageTk.PhotoImage(image = self.slice_pil_Enface)
    self.image_on_canvas_Enface = self.canvasXYCorrection.create_image(0, 0, image=self.photoEnface, anchor=NW)        
    
    self.vol_sliderXY_corr.grid(row=0,column=0, padx=pad_x_sliderXY,pady=(1,0),sticky="NW")
    self.vol_sliderXY_corr.set(self.volume.shape[0]//2)
    if self.rectid == None:
        _createVariables(self.master)
        _createCanvasBinding()
    
    #check if GT exists
    if self.GT_EXISTS:
        createHeatmap()
    
    self.canvasXZCorrection.grid(row=1,column=0, columnspan = 3, padx=(200,0), pady=(20,0))
    #mouse bindings for correction
    self.canvasXZCorrection.bind("<Button-1>", leftclickXYCorr)
    self.canvasXZCorrection.bind('<B1-Motion>', leftclickXYCorrMovement) 
    
    self.sbarV_XZCorr = Scrollbar(self.correctionFrame, orient=VERTICAL,command=self.canvasXZCorrection.yview)
    self.canvasXZCorrection.config(yscrollcommand=self.sbarV_XZCorr.set)
    self.sbarV_XZCorr.grid(row=1, column = 3, sticky=N+S,pady=(20,0))
    self.canvasXZCorrection.configure(scrollregion = self.canvasXZCorrection.bbox("all"))

    self.sliceXZCorrection = self.volumeXZ[self.volumeXZ.shape[0]//2].astype('uint8')
    self.sliceXZCorrection  = rsc(self.sliceXZCorrection , 3, mode='reflect', preserve_range = True).astype('uint8')
    self.slice_pilXZCorrection = PIL.Image.fromarray(self.sliceXZCorrection)
    self.photoXZCorrection = PIL.ImageTk.PhotoImage(image = self.slice_pilXZCorrection)
    self.image_on_canvas_XZCorrection = self.canvasXZCorrection.create_image(0,0, image=self.photoXZCorrection, anchor=NW)
    
    self.vol_sliderXZCorrection = Scale(self.correctionFrame, from_=0, length= 350, to=self.volumeXZ.shape[0]-1, orient=VERTICAL, command=updateVolumeXZCorrSlice,bg=self.canvasbackground, fg =self.canvasforeground)
    self.vol_sliderXZCorrection.grid(row=1,column=4,sticky=W,pady=(20,0))
    self.vol_sliderXZCorrection.set(self.volumeXZ.shape[0]//2)  

    ctr_row = 0
    self.statusText = StringVar()
    self.statusText.set("")
    self.LabelPGBarCorr = Label(self.correctionFrame, width = 18,height= 2, textvariable=self.statusText, font=self.clrbarFont,anchor="nw", bg=self.canvasbackground, fg =self.hlghtbg)
    self.LabelPGBarCorr.grid(row=0,column=0,pady=(10,0),padx=(5,0),sticky='NW')
    
    #button menu - most callbacks is put in threading to not stop the GUI
    self.btnExplore = Button(self.correctionFrame,width = 15 , height = self.manCorrBtnHeight, text="Explore Mode",font=self.statusFont, command=self.exploremode, activebackground=self.btn_common_bg, fg=self.btnforeground, bg=self.btnbackground, activeforeground=self.btnforeground)
    self.btnExplore.grid(row=ctr_row,column=0,padx=(40,0),pady=(self.corr_y_pad_increment,0), sticky="NW")
    CreateToolTip(self.btnExplore, self.ttip_dict['manExplore'])
    self.btnLoadGT = Button(self.correctionFrame,width = 15 , height = self.manCorrBtnHeight, text="Load Ground Truth",font=self.statusFont, command=lambda:background(loadGroundtruth), activebackground=self.btn_common_bg, fg=self.btnforeground, bg=self.btnbackground, activeforeground=self.btnforeground)
    self.btnLoadGT.grid(row=ctr_row,column=0,padx=(40,0),pady=(self.corr_y_pad_increment*2,0), sticky="NW")
    CreateToolTip(self.btnLoadGT, self.ttip_dict['manLoadGT'])
    self.btnSaveRes = Button(self.correctionFrame,width = 15 , height = self.manCorrBtnHeight, text="Save Result",font=self.statusFont, command=lambda:background(saveSegmentation), activebackground=self.btn_common_bg, fg=self.btnforeground, bg=self.btnbackground, activeforeground=self.btnforeground)
    self.btnSaveRes.grid(row=ctr_row,column=0,padx=(40,0),pady=(self.corr_y_pad_increment*3,0), sticky="NW")
    CreateToolTip(self.btnSaveRes, self.ttip_dict['manSaveRes'])
    '''
        Drop Down Menu to select correction layer
        IntegrateOption: If new layers are added, integrate it here by adding the specific layer.
    '''
    Label(self.correctionFrame, width = 11,height= 1, text="Layer",font=self.headerFont,anchor="nw", bg=self.canvasbackground, fg =self.canvasforeground).grid(row=ctr_row,column=0,padx=(40,0),pady=(self.corr_y_pad_increment*4,0),sticky="NW")
    self.optionManCorrLayer = OptionMenu(self.correctionFrame, self.varManCorrLayer, "ILM", "RPE", "BM", command=layerChange)
    self.optionManCorrLayer.config(font=self.textFont)
    menu = self.optionManCorrLayer.nametowidget(self.optionManCorrLayer.menuname)
    menu.configure(font=self.textFont)
    self.optionManCorrLayer.grid(row=ctr_row,column=0,padx=(40,0),pady=(int(self.corr_y_pad_increment*4.5),0),sticky="NW")
    CreateToolTip(self.optionManCorrLayer, self.ttip_dict['manCorrLayer'])

    ctr_row =+1
    self.btnRunManRefine = Button(self.correctionFrame,width = 15 , height = self.manCorrBtnHeight + 1, text="Run Refinement",font=self.statusFont, command=lambda:background(propagateCorrection), activebackground=self.btn_common_bg, fg=self.btnforeground, bg=self.btnbackground, activeforeground=self.btnforeground)
    self.btnRunManRefine.grid(row=ctr_row,column=0,padx=(40,0),pady=(22,0), sticky="NW")
    CreateToolTip(self.btnRunManRefine, self.ttip_dict['manRunRefine'])
    self.btnResetSlice = Button(self.correctionFrame,width = 15 , height = self.manCorrBtnHeight + 1, text="Reset slice",font=self.statusFont, command=lambda:background(resetSlice), activebackground=self.btn_common_bg, fg=self.btnforeground, bg=self.btnbackground, activeforeground=self.btnforeground)
    self.btnResetSlice.grid(row=ctr_row,column=0,padx=(40,0),pady=(305,0), sticky="NW")
    CreateToolTip(self.btnResetSlice, self.ttip_dict['manResetSlice'])

    self.canvasXZCorrection.yview_moveto(0.5)