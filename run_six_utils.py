import sys
sys.path.append('/home/s/seanmacb/Camera/pkgs')
sys.path.append('/home/s/seanmacb/Camera/pkgs/eo_pipe/python/lsst/eo/pipe')
import myutils_ar.myutils as myu
import eo_pipe.python.lsst.eo.pipe as eo_pipe
from eo_pipe.python.lsst.eo.pipe import (readNoiseTask, darkCurrentTask, defectsTask, eperTask, divisaderoTearingTask, ptcPlotsTask,linearityPlotsTask, bfAnalysisTask)
import lsst.daf.butler as daf_butler
from lsst.obs.lsst import LsstCam, LsstTS8
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1 import make_axes_locatable
import lsst.geom
from lsst.afw import cameraGeom
import lsst.afw.display as afwDisplay
afwDisplay.setDefaultBackend('matplotlib') 
import pandas as pd
import matplotlib.ticker as ticker
from scipy.stats import norm
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from lsst.ip.isr.isrTask import IsrTask
import lsst.afw.image as afwImage
from lsst.afw.image import Mask
import lsst.geom as geom
from matplotlib.lines import Line2D

### START NEW STUFF

class defectMaskComparison:
    def __init__(self, detector_num,collections='LSSTCam/raw/all,LSSTCam/calib'.split(",")):
        self.detector_num = detector_num
        self.BBox = getDetectorBBox([self.detector_num])[0] # remove rsu here
        # self.detector = detector_for_Defects
        self.collections = collections
        self.runList = []
        return
    
    def setRun(self,run,collection = None):
        self.currentRun = run
        if type(self.currentRun)==str:
            self.repo="embargo_new"
        elif self.currentRun>13999:
            self.repo = "/repo/main"     
        else:
            self.repo = "/repo/ir2"

        if collection != None:
            self.setCollections(collection)
    
        self.butler = makeButler(self.repo,self.collections) # remove rsu here
        self.registry = self.butler.registry
        self.runList.append(self.currentRun)
        return

    def setCollections(self,collections):
        self.collections = collections
        return
    def makeMaskedDetectorImage(self):
        self.maskedDetectorImage = afwImage.MaskedImageF(self.BBox)
        return
    def setDetectorObj(self,detector_bay_slot):
        self.bay_slot = detector_bay_slot
        self.detectorObject = myu.camera[self.bay_slot]
        return
    def setAmplifier(self,segment):
        self.segment = segment
        self.amplifierObject = self.detectorObject[self.segment]
        return
    def setAmplifierBBox(self):
        self.ampBBox = self.amplifierObject.getBBox()
        return
    def combineMasks(self,masks,runList):
        # masks is an array of defectMaskComparison objects
        
        # Scale masks properly
        self.planeDefs = {"None":0}
        for mask,num,maskName in zip(masks,range(len(masks)),runList):
            mask *= 2**num
            self.planeDefs[str(maskName)] = 2**num 
        # Add masks together
        self.combinedMask = masks.pop()
        for mask in masks:
            self.combinedMask+=mask
        # Need to add something here to add additional key, value pairs to the self.planeDefs dict
        self.setMaskX()
        self.makeMaskedAmpImage()
        return
        
    def setMaskX(self,adjustment=0):
        minx,maxx,miny,maxy = self.amplifierObject.getBBox().getMinX(),self.amplifierObject.getBBox().getMaxX(),self.amplifierObject.getBBox().getMinY(),self.amplifierObject.getBBox().getMaxY() # Changed mask.amplifierObject to self.amplifierObject
        self.ampMaskArray = self.combinedMask[miny:maxy+1,minx:maxx+1]
        self.updateplaneDefs(adjustment=adjustment)
        self.ampMaskX = Mask(geom.Box2I(self.amplifierObject.getBBox()),self.planeDefs) # changed mask.planeDefs to self.planeDefs
        self.ampMaskX.array = self.ampMaskArray
        
        self.ampMaskXForMTV = Mask(geom.Box2I(self.amplifierObject.getBBox()),self.planeDefs) # changed mask.planeDefs to self.planeDefs
        self.ampMaskXForMTV.array = self.ampMaskArray
        self.ampMaskXForMTV.array[np.where(self.ampMaskXForMTV.array!=0)] = 2**self.ampMaskXForMTV.array[np.where(self.ampMaskXForMTV.array!=0)]
        
        # minx,maxx,miny,maxy = self.amplifierObject.getBBox().getMinX(),self.amplifierObject.getBBox().getMaxX(),self.amplifierObject.getBBox().getMinY(),self.amplifierObject.getBBox().getMaxY()
        # self.ampMaskX = Mask(self.combinedMask[miny:maxy+1,minx:maxx+1],planeDefs=self.planeDefs)

    def updateplaneDefs(self,adjustment=0):
        newMaskDefs,newMaskDefs_filtered = {},{}
        for num in range(2**len(self.runList)):
            bin_num = format(num, "0{}b".format(len(self.runList)))
            first=False
            binary_list_string = list(bin_num)
            binary_list_ints = [int(s) for s in binary_list_string]
            
            key = ""
            # print("binary_list: {}".format(binary_list_ints))
            for bin_number,iterator in zip(binary_list_ints,range(len(binary_list_ints)-1,-1,-1)): # This is where I updated the plane definitions
                if iterator!=len(binary_list_ints)-1 and key!="" and bin_number==1: # I changed the iterator condition from iterator!=0 to iterator!=len(binary_list_ints)-1
                    key+=", "
                if bin_number==1:
                    key+=str(self.runList[iterator])
            if key == "":
                key = "None"
    ###
            # if key!="None":
            newMaskDefs[key] = num + adjustment
    ###
            if num in self.ampMaskArray:
                newMaskDefs_filtered[key] = num + adjustment

        self.minimalPlaneDefs = newMaskDefs_filtered
        self.planeDefs = newMaskDefs
        

    def makeMaskedAmpImage(self):
        self.ampImage = afwImage.MaskedImageF(self.amplifierObject.getBBox(),self.planeDefs)
        self.ampImage.setMask(self.ampMaskXForMTV)

        
    def makeMask(self,maskColl):
        if maskColl==None:
            self.entry_collection = self.registry.queryCollections(f"*{self.defectPath}*{self.currentRun}*",collectionTypes=daf_butler.CollectionType.CHAINED)
        else:
            self.entry_collection = maskColl
        # self.entry_collection = self.registry.queryCollections(f"*{self.defectPath}*{self.currentRun}*",collectionTypes=daf_butler.CollectionType.CHAINED)
        self.maskKwargs = makeKwargs([self.currentRun],self.defectDatasetType,self.detector_num,self.entry_collection) # remove rsu here
        # I implemented this to catch oddities in the makeKwargs function, but those oddities are actually just that eo_pipe was not ran for run 5 EO runs
        # if len(mask3.maskKwargs['collections'])==0: 
        #     mask3.maskKwargs['collections'] = self.collections
        self.maskRefs = list(self.registry.queryDatasets(**self.maskKwargs)) # here is where it is not doing the proper thing...

        found=False
        for ref in self.maskRefs:
            if self.entry_collection+"/" in ref.run:
                self.maskedDefect = self.butler.get(ref)
                print("For run {}, ref id: {}".format(self.currentRun,ref.dataId))
                found=True
                break
        if not found:
            print("Ref not found for {}".format(self.maskKwargs))
            print("Repo: {}".format(self.repo))
        # return -1 # If not found
        
        # if len(self.maskRefs)==1:
        #     self.maskedDefect = self.butler.get(self.maskRefs[0])

    def selectMaskedDefect(self,ref):
        self.maskedDefect = self.butler.get(ref)
        
    def applyMask(self,image):
        self.maskedDefect.maskPixels(image) # add maskName=datasetType here

    def setDefectPath(self,defectPath):
        self.defectPath = defectPath

    def setDefectDatasetType(self,defectDatasetType):
        self.defectDatasetType = defectDatasetType

    def applyMaskToImage(self,defectPath,defectDatasetType,maskColl=None):
        ## Use this function, which calls the above
        self.setDefectPath(defectPath)
        self.setDefectDatasetType(defectDatasetType)
        self.makeMask(maskColl=maskColl)
        if "maskedDefect" not in dir(self):
            for ref,num in zip(self.maskRefs,range(len(self.maskRefs))):
                print(num,ref.run)
            iterator = int(input("Please select a ref from the ref-set\n"))
            self.selectMaskedDefect(self.maskRefs[iterator])
            
        self.makeMaskedDetectorImage()
        self.applyMask(self.maskedDetectorImage)
        self.maskArray = self.maskedDetectorImage.mask.array

    def adjustMaskValues(self,adjustment):
        self.combinedMask[np.where(self.combinedMask!=0)]+=adjustment
        self.setMaskX(adjustment=adjustment)
        self.makeMaskedAmpImage()

    def scaleMaskForMTV(self):
        self.ampMaskX.array[np.where(self.ampMaskX.array!=0)] = 2**self.ampMaskX.array[np.where(self.ampMaskX.array!=0)]

def showAndClear():
    plt.show()
    # Clear the current axes.
    plt.cla() 
    # Clear the current figure.
    plt.clf() 
    # Closes all the figure windows.
    plt.close('all')   
    plt.close(fig)
    
    return

def oppoImage(defectDatasetType):
    if defectDatasetType.__contains__("right"):
        return "dark"
    else:
        return "flat"

def defectStabilityFig(mask,defectDatasetType,yLabel,xlim=None,ylim=None,save=False,figure_dir=None,clim=None):
    fig,axs = plt.subplots(1,3,figsize=[16,8])

    fig.suptitle("{}_{}\n{}".format(mask.bay_slot,mask.segment,defectDatasetType))

    # Image of the region, with colorbar
    myImgType = oppoImage(defectDatasetType)
    # Get the image here
    img = mask.butler.get(list(mask.butler.registry.queryDatasets(myImgType,where="detector={}".format(mask.detector_num)))[0])

    for ax in axs[0:2]:
        if ylim == None and xlim ==None:
            img = img[mask.amplifierObject.getBBox()]
        else: 
            if ylim!=None:
                ax.set_ylim(ylim)
                img = img[:,ylim[0]:ylim[1]]
                mask.ampImage = mask.ampImage[:,ylim[0]:ylim[1]]
            else:
                ##
                img = img[:,mask.amplifierObject.getBBox().minY:mask.amplifierObject.getBBox().maxY]
            if xlim !=None:
                ax.set_xlim(xlim)
                img = img[xlim[0]:xlim[1],:]
                mask.ampImage = mask.ampImage[xlim[0]:xlim[1],:]
            else:
                ##
                img = img[mask.amplifierObject.getBBox().minX:mask.amplifierObject.getBBox().maxX,:]
                
    plt.sca(axs[0])
    disp = afwDisplay.Display(fig)
    disp.setImageColormap('cool')
    # mycolors = format_mtv(disp,mask,axs[0])
    myMin,myMax = np.nanpercentile(img.image.array,[10,90])
    disp.scale('linear',myMin,myMax)
    disp.mtv(img)
    axs[0].set_title("{}, run {}".format(myImgType,mask.currentRun))
    
# Subplot of masked region of detector, on amp level

    plt.sca(axs[1])
    myDisp = afwDisplay.Display(fig)
    myDisp.setImageColormap('bwr')
    mycolors = format_mtv(myDisp,mask,axs[1])
    myDisp.mtv(mask.ampImage)
    leg_elements,mycolors2 = formatmtv2(myDisp,mask,axs[0],mycolors)
    axs[1].legend(handles=leg_elements,loc=(1.1,0.9))
    myDisp.show_colorbar(False)
    # disp.scale('linear',0,1)
        
    overlaps,n_overlaps = eo_stats(mask,amplifier=True)
    
    # Subplot of the changes over multiple different EO runs   
    runColorDict = getColorsRunPlot(disp,mask)
    plt.sca(axs[2])
    axs[2].scatter(np.arange(len(overlaps)),n_overlaps,color=mycolors2)
    format_runPlot(axs[2],leg_elements,yLabel)

    # ax.set_title()
    fig.tight_layout()
    if save:
        plt.savefig(figure_dir+"/{}_{}.jpg".format(mask.bay_slot,mask.segment),dpi=200)

    showAndClear(fig)
    del fig
    del myDisp
    del img
    return 

def maskMerger(detector,seg,runList,repo_coll,defectPath,defectDatasetType,maskColls=None):
    # Do setup for this amp
    myMaskArrayList = []
    myMaskList = []
    myMask = defectMaskComparison(myu.det_nums[detector])
    myMask.setDetectorObj(detector)
    myMask.setAmplifier(seg)
    for run,repo_collec,maskColl in zip(runList,repo_coll,maskColls):
        repo = repo_collec[0]
        collection = repo_collec[1]
        
        myMask.setRun(run,collection) # This function looks fine
        
        myMask.applyMaskToImage(defectPath,defectDatasetType,maskColl=maskColl)
    
        myMaskArrayList.append(myMask.maskArray)
        
    myMask.combineMasks(myMaskArrayList,myMask.runList)
    return myMask

def eo_stats(mask,amplifier=False,verbose=False):
    All_combos = np.unique(mask.ampMaskX.array)
    
    flattened_mask = np.reshape(mask.ampMaskX.array,(1,-1))[0]
    if verbose:
        print("Detector {}, {}, {}".format(mask.detector_num,myu.det_names[mask.detector_num],get_serial_from_number([mask.detector_num])[0])) # remove rsu from get_serial_from_number
        if amplifier:
            print("Segment {}".format(mask.segment))
            print("Total pixels in amplifier: {}".format(len(flattened_mask)))
    
    overlaps = []
    n_overlaps = []

    for number in np.arange(2**len(mask.runList)+1):
        binary = format(number, "0{}b".format(len(mask.runList)))
        binary_list_string = list(binary)
        binary_list_ints = [int(s) for s in binary_list_string]
        runListOverlaps = []
        # for bin_number,iterator in zip(binary_list_ints,range(len(binary_list_ints)-1,-1,-1)):
        for bin_number,iterator in zip(binary_list_ints,range(len(binary_list_ints))):
            if bin_number==1:
                runListOverlaps.append(mask.runList[iterator])
        number_overlaps = len(np.where(flattened_mask == number)[0])
        
        if number in All_combos:
            if number==0:
                base_str = "Number of pixels that are not masked in any EO runs"
            else:
                base_str = "Number of pixels that overlap from EO runs {}".format(runListOverlaps)
                overlaps.append(runListOverlaps)
                n_overlaps.append(number_overlaps)
            
        if verbose:        
            print("{}: {}".format(base_str,number_overlaps))
            
        if verbose:
            print("Number of pixels that are masked in any EO run: {}".format(len(np.where(flattened_mask > 0)[0])),end="\n\n")
    
    return overlaps,n_overlaps

colorList = ["#ff0000","#00ff00","#0000ff","#ffa0f7","#918BE1","#fea60e",
             "#161f7e","#6600cc","#ff00ff","#ffff00","#00ffff","#66ffff"]

def format_mtv(disp,mask,ax):
    base_dic = {'BAD': afwDisplay.IGNORE,'CR': afwDisplay.IGNORE,'EDGE': afwDisplay.IGNORE,'INTERPOLATED': afwDisplay.IGNORE,'SATURATED': afwDisplay.IGNORE,'DETECTED': afwDisplay.IGNORE,'DETECTED_NEGATIVE': afwDisplay.IGNORE,'SUSPECT': afwDisplay.IGNORE,'NO_DATA': afwDisplay.IGNORE,'INTRP': afwDisplay.IGNORE,'SAT': afwDisplay.IGNORE}
    disp.setMaskPlaneColor(base_dic)
    mycolors = []
    color_iter = 0
    for key in list(mask.minimalPlaneDefs.keys()):
        if key!="None":
            # print(key,colorList[color_iter])
            disp.setMaskPlaneColor(key,colorList[color_iter])
            disp.setMaskTransparency(0,name=key)
            mycolors.append(colorList[color_iter])
            color_iter+=1
        # else:
        #     print(key,afwDisplay.IGNORE)
        #     disp.setMaskPlaneColor(key,afwDisplay.IGNORE)
        #     disp.setMaskTransparency(0,name=key)
    ax.set_xlabel("Serial Register")
    ax.set_ylabel("Parallel Register")


    return mycolors
        # else:
        #     disp.setMaskPlaneColor(key,'black')
        # else:
            
    
    # startingDict = disp.getMaskPlaneColor()
    
    # if len(startingDict)<len(mask.planeDefs):
    #     for color,key in zip(colorList,len(mask.planeDefs) - len(startingDict)):
    #         disp.setMaskPlaneColor(key,color)
    #         disp.setMaskTransparency(0,key)
    
def formatmtv2(disp,mask,ax,mycolors):
    legend_elements = []

    startingDict = disp.getMaskPlaneColor()
    mycolors2 = []
    color_iter=0
    for key,mask_key in zip(list(startingDict.keys()),list(mask.planeDefs.keys())):
        if mask_key in list(mask.minimalPlaneDefs.keys()):
            # disp.setMaskPlaneColor(key,colorList[color_iter])
            # disp.setMaskTransparency(100,key)
            if mask_key!="None":
                legend_elements.append(Line2D([0], [0], marker='s', color=mycolors[color_iter], label=mask_key,
                                                  markerfacecolor=mycolors[color_iter], markersize=8,ls="None"))
                mycolors2.append(mycolors[color_iter])
                color_iter+=1

    return legend_elements,mycolors2

def getColorsRunPlot(disp,mask):
    starterDict = {}
    iter_ = 0
    for key,value in zip(list(disp.getMaskPlaneColor().keys()),list(disp.getMaskPlaneColor().values())):
        if key in mask.minimalPlaneDefs and key!="None":
            starterDict[key] = colorList[iter_]
            iter_ +=1
    return starterDict

def format_runPlot(ax,overlaps,ylabel):
    ax.grid()
    if len(overlaps)!=0:
        ax.set_xticks(np.arange(len(overlaps)))
        ax.set_xlim(np.min(np.arange(len(overlaps)))-0.5,np.max(np.arange(len(overlaps)))+0.5)
        ax.xaxis.set_ticklabels([str(s).split("(")[1].split(")")[0] for s in overlaps], rotation=45, fontsize=10)
    ax.set_ylabel(ylabel)
    return

###

### OLD STUFF

# class defectMaskComparison:
#     def __init__(self, detector_num,collections='LSSTCam/raw/all,LSSTCam/calib'.split(",")):
#         self.detector_num = detector_num
#         self.BBox = getDetectorBBox([self.detector_num])[0]
#         # self.detector = detector_for_Defects
#         self.collections = collections
#         self.runList = []
    
#     def setRun(self,run):
#         self.currentRun = run
#         if type(self.currentRun)==np.str_ or type(self.currentRun)==str:
#             self.repo="embargo_new"
#         elif self.currentRun>13999:
#             self.repo = "/repo/main"     
#         else:
#             self.repo = "/repo/ir2"
#         self.butler = makeButler(self.repo,self.collections)
#         self.registry = self.butler.registry
#         self.runList.append(self.currentRun)
        
#     def makeMaskedDetectorImage(self):
#         self.maskedDetectorImage = afwImage.MaskedImageF(self.BBox)

#     def setDetectorObj(self,detector_bay_slot):
#         self.bay_slot = detector_bay_slot
#         self.detectorObject = myu.camera[self.bay_slot]

#     def setAmplifier(self,segment):
#         self.segment = segment
#         self.amplifierObject = self.detectorObject[self.segment]

#     def setAmplifierBBox(self):
#         self.ampBBox = self.amplifierObject.getBBox()
    
#     def combineMasks(self,masks,runList):
#         # masks is an array of defectMaskComparison objects
        
#         # Scale masks properly
#         self.planeDefs = {"None":0}
#         for mask,num,maskName in zip(masks,range(len(masks)),runList):
#             mask *= 2**num
#             self.planeDefs[str(maskName)] = 2**num 
#         # Add masks together
#         self.combinedMask = masks.pop()
#         for mask in masks:
#             self.combinedMask+=mask
#         # Need to add something here to add additional key, value pairs to the self.planeDefs dict
#         self.setMaskX()
#         self.makeMaskedAmpImage()

#     def setMaskX(self,adjustment=0):
#         minx,maxx,miny,maxy = self.amplifierObject.getBBox().getMinX(),self.amplifierObject.getBBox().getMaxX(),self.amplifierObject.getBBox().getMinY(),self.amplifierObject.getBBox().getMaxY()
#         self.ampMaskArray = self.combinedMask[miny:maxy+1,minx:maxx+1]
#         self.updateplaneDefs(adjustment=adjustment)
#         self.ampMaskX = Mask(geom.Box2I(self.amplifierObject.getBBox()),self.planeDefs)
#         self.ampMaskX.array = self.ampMaskArray
        
#         self.ampMaskXForMTV = Mask(geom.Box2I(self.amplifierObject.getBBox()),self.planeDefs)
#         self.ampMaskXForMTV.array = self.ampMaskArray
#         self.ampMaskXForMTV.array[np.where(self.ampMaskXForMTV.array!=0)] = 2**self.ampMaskXForMTV.array[np.where(self.ampMaskXForMTV.array!=0)]
        
#         # minx,maxx,miny,maxy = self.amplifierObject.getBBox().getMinX(),self.amplifierObject.getBBox().getMaxX(),self.amplifierObject.getBBox().getMinY(),self.amplifierObject.getBBox().getMaxY()
#         # self.ampMaskX = Mask(self.combinedMask[miny:maxy+1,minx:maxx+1],planeDefs=self.planeDefs)

#     def updateplaneDefs(self,adjustment=0):
#         newMaskDefs,newMaskDefs_filtered = {},{}
#         for num in range(2**len(self.runList)):
#             bin_num = format(num, "0{}b".format(len(self.runList)))
#             first=False
#             binary_list_string = list(bin_num)
#             binary_list_ints = [int(s) for s in binary_list_string]
            
#             key = ""
#             for bin_number,iterator in zip(binary_list_ints,range(len(binary_list_ints))):
#                 if iterator!=0 and key!="" and bin_number==1:
#                     key+=", "
#                 if bin_number==1:
#                     key+=str(self.runList[iterator])
#             if key == "":
#                 key = "None"
#     ###
#             # if key!="None":
#             newMaskDefs[key] = num + adjustment
#     ###
#             if num in self.ampMaskArray:
#                 newMaskDefs_filtered[key] = num + adjustment

#         self.minimalPlaneDefs = newMaskDefs_filtered
#         self.planeDefs = newMaskDefs
        

#     def makeMaskedAmpImage(self):
#         self.ampImage = afwImage.MaskedImageF(self.amplifierObject.getBBox(),self.planeDefs)
#         self.ampImage.setMask(self.ampMaskXForMTV)

        
#     def makeMask(self,maskColl=None):
#         if maskColl==None:
#             self.entry_collection = self.registry.queryCollections(f"*{self.defectPath}*{self.currentRun}*",collectionTypes=daf_butler.CollectionType.CHAINED)
#         else:
#             self.entry_collection = maskColl
#         self.maskKwargs = makeKwargs([self.currentRun],self.defectDatasetType,self.detector_num,self.entry_collection)
#         # I implemented this to catch oddities in the makeKwargs function, but those oddities are actually just that eo_pipe was not ran for run 5 EO runs
#         # if len(mask3.maskKwargs['collections'])==0: 
#         #     mask3.maskKwargs['collections'] = self.collections
#         self.maskRefs = list(self.registry.queryDatasets(**self.maskKwargs))
#         if len(self.maskRefs)==1:
#             self.maskedDefect = self.butler.get(self.maskRefs[0])

#     def selectMaskedDefect(self,ref):
#         self.maskedDefect = self.butler.get(ref)
        
#     def applyMask(self,image):
#         self.maskedDefect.maskPixels(image) # add maskName=datasetType here

#     def setDefectPath(self,defectPath):
#         self.defectPath = defectPath

#     def setDefectDatasetType(self,defectDatasetType):
#         self.defectDatasetType = defectDatasetType

#     def applyMaskToImage(self,defectPath,defectDatasetType,maskColl=None,override=None):
#         ## Use this function, which calls the above
#         self.setDefectPath(defectPath)
#         self.setDefectDatasetType(defectDatasetType)
#         self.makeMask(maskColl=maskColl)
#         if "maskedDefect" not in dir(self):
#             if override==None:
#                 for ref,num in zip(self.maskRefs,range(len(self.maskRefs))):
#                     print(num,ref.run)
#                 iterator = int(input("Please select a ref from the ref-set\n"))
#                 self.selectMaskedDefect(self.maskRefs[iterator])
#             else:
#                 self.selectMaskedDefect(self.maskRefs[override])
            
#         self.makeMaskedDetectorImage()
#         self.applyMask(self.maskedDetectorImage)
#         self.maskArray = self.maskedDetectorImage.mask.array

#     def adjustMaskValues(self,adjustment):
#         self.combinedMask[np.where(self.combinedMask!=0)]+=adjustment
#         self.setMaskX(adjustment=adjustment)
#         self.makeMaskedAmpImage()

#     def scaleMaskForMTV(self):
#         self.ampMaskX.array[np.where(self.ampMaskX.array!=0)] = 2**self.ampMaskX.array[np.where(self.ampMaskX.array!=0)]

# def showAndClear():
#     plt.show()
#     # Clear the current axes.
#     plt.cla() 
#     # Clear the current figure.
#     plt.clf() 
#     # Closes all the figure windows.
#     plt.close('all')   
#     plt.close(fig)
    
#     return

# def eo_stats(mask,amplifier=False,verbose=False):
#     All_combos = np.unique(mask.ampMaskX.array)
    
#     flattened_mask = np.reshape(mask.ampMaskX.array,(1,-1))[0]
#     if verbose:
#         print("Detector {}, {}, {}".format(mask.detector_num,myu.det_names[mask.detector_num],get_serial_from_number([mask.detector_num])[0]))
#         if amplifier:
#             print("Segment {}".format(mask.segment))
#             print("Total pixels in amplifier: {}".format(len(flattened_mask)))
    
#     overlaps = []
#     n_overlaps = []

#     for number in np.arange(2**len(mask.runList)+1):
#         binary = format(number, "0{}b".format(len(mask.runList)))
#         binary_list_string = list(binary)
#         binary_list_ints = [int(s) for s in binary_list_string]
#         runListOverlaps = []
#         for bin_number,iterator in zip(binary_list_ints,range(len(binary_list_ints))):
#             if bin_number==1:
#                 runListOverlaps.append(mask.runList[iterator])
#         number_overlaps = len(np.where(flattened_mask == number)[0])
        
#         if number in All_combos:
#             if number==0:
#                 base_str = "Number of pixels that are not masked in any EO runs"
#             else:
#                 base_str = "Number of pixels that overlap from EO runs {}".format(runListOverlaps)
#                 overlaps.append(runListOverlaps)
#                 n_overlaps.append(number_overlaps)
            
#         if verbose:        
#             print("{}: {}".format(base_str,number_overlaps))
            
#         if verbose:
#             print("Number of pixels that are masked in any EO run: {}".format(len(np.where(flattened_mask > 0)[0])),end="\n\n")
    
#     return overlaps,n_overlaps

# colorList = ["#ff0000","#00ff00","#0000ff","#ffa0f7","#918BE1","#fea60e",
#              "#161f7e","#6600cc","#ff00ff","#ffff00","#00ffff","#66ffff"]

# def format_mtv(disp,mask,ax):
#     base_dic = {'BAD': afwDisplay.IGNORE,'CR': afwDisplay.IGNORE,'EDGE': afwDisplay.IGNORE,'INTERPOLATED': afwDisplay.IGNORE,'SATURATED': afwDisplay.IGNORE,'DETECTED': afwDisplay.IGNORE,'DETECTED_NEGATIVE': afwDisplay.IGNORE,'SUSPECT': afwDisplay.IGNORE,'NO_DATA': afwDisplay.IGNORE,'INTRP': afwDisplay.IGNORE,'SAT': afwDisplay.IGNORE}
#     disp.setMaskPlaneColor(base_dic)
#     mycolors = []
#     color_iter = 0
#     for key in list(mask.minimalPlaneDefs.keys()):
#         if key!="None":
#             # print(key,colorList[color_iter])
#             disp.setMaskPlaneColor(key,colorList[color_iter])
#             disp.setMaskTransparency(0,name=key)
#             mycolors.append(colorList[color_iter])
#             color_iter+=1
#         # else:
#         #     print(key,afwDisplay.IGNORE)
#         #     disp.setMaskPlaneColor(key,afwDisplay.IGNORE)
#         #     disp.setMaskTransparency(0,name=key)
#     ax.set_xlabel("Serial Register")
#     ax.set_ylabel("Parallel Register")


#     return mycolors
#         # else:
#         #     disp.setMaskPlaneColor(key,'black')
#         # else:
            
    
#     # startingDict = disp.getMaskPlaneColor()
    
#     # if len(startingDict)<len(mask.planeDefs):
#     #     for color,key in zip(colorList,len(mask.planeDefs) - len(startingDict)):
#     #         disp.setMaskPlaneColor(key,color)
#     #         disp.setMaskTransparency(0,key)
    
# def formatmtv2(disp,mask,ax,mycolors):
#     legend_elements = []

#     startingDict = disp.getMaskPlaneColor()
#     mycolors2 = []
#     color_iter=0
#     for key,mask_key in zip(list(startingDict.keys()),list(mask.planeDefs.keys())):
#         if mask_key in list(mask.minimalPlaneDefs.keys()):
#             # disp.setMaskPlaneColor(key,colorList[color_iter])
#             # disp.setMaskTransparency(100,key)
#             if mask_key!="None":
#                 legend_elements.append(Line2D([0], [0], marker='s', color=mycolors[color_iter], label=mask_key,
#                                                   markerfacecolor=mycolors[color_iter], markersize=8,ls="None"))
#                 mycolors2.append(mycolors[color_iter])
#                 color_iter+=1

#     return legend_elements,mycolors2

# def getColorsRunPlot(disp,mask):
#     starterDict = {}
#     iter_ = 0
#     for key,value in zip(list(disp.getMaskPlaneColor().keys()),list(disp.getMaskPlaneColor().values())):
#         if key in mask.minimalPlaneDefs and key!="None":
#             starterDict[key] = colorList[iter_]
#             iter_ +=1
#     return starterDict

# def format_runPlot(ax,overlaps,ylabel):
#     ax.grid()
#     ax.set_xticks(np.arange(len(overlaps)))
#     ax.set_xlim(np.min(np.arange(len(overlaps)))-0.5,np.max(np.arange(len(overlaps)))+0.5)
#     ax.xaxis.set_ticklabels([str(s)[1:-1] for s in overlaps], rotation=45, fontsize=10)
#     ax.set_ylabel(ylabel)
#     return

### END OLD STUFF

camera = myu.LsstCam.getCamera()

def showAndClear(fig):
    plt.show()
    # Clear the current axes.
    plt.cla() 
    # Clear the current figure.
    plt.clf() 
    # Closes all the figure windows.
    plt.close('all')   
    plt.close(fig)
    plt.close()
    return

def MakeHistogram(img,
                  run,
                  detector,
                  run_num,
                  face_color='#2ab0ff',
                  bin_num=100,
                  edgecolor='#e0e0e0',
                  lw=0.5,
                  alpha=0.7,
                  xlim=(0,5E4),
                  y_ax='linear',
                  ylim=(1E-10,1E-3),
                  fit_norm=True,
                  n_sig=1):
    image = getImage(img,run,detector)
    
    arr = image.image.array.flatten()

    fig,ax = plt.subplots(figsize=[8,6])

    if y_ax == 'log':
        ax.semilogy()
    elif y_ax != 'linear':
        print("y_ax keyword not set correctly")
        return -1        
    
    n, bins, patches = ax.hist(arr, bins=bin_num, facecolor=face_color, edgecolor=edgecolor, linewidth=lw, alpha=alpha,density=True)

    n = n.astype('int') # it MUST be integer# Good old loop. Choose colormap of your taste
    
    # for i in range(len(patches)):
    #     patches[i].set_facecolor(plt.cm.plasma(n[i]/max(n)))# Make one bin stand out   
    
    # patches[47].set_fc('red') # Set color
    # patches[47].set_alpha(1) # Set opacity# Add annotation
    # plt.annotate('Important Bar!', xy=(0.57, 175), xytext=(2, 130), fontsize=15, arrowprops={'width':0.4,'headwidth':7,'color':'#333333'})# Add title and labels with custom font sizes
    
    ax.set_title('Distribution of pixel values from detector {var1} ({var3}) {var4} image, run {var2}'.format(var1=detector,var2=run_num,var3=get_serial_from_number([detector])[0],var4=img), fontsize=12)
    ax.set_xlabel('Counts [e-]', fontsize=10)
    ax.set_ylabel('Normalized counts [# of pixels]', fontsize=10)


    if fit_norm:
        mu_kwargs = {
            "color":"green",
            "ls":"--"
        }
        sigma_kwargs = {
            'color':"blue"
        }
        # Fit a normal distribution to the data:
        mu, std = norm.fit(arr)
        # Plot the PDF.
        xmin, xmax = xlim
        x = np.linspace(xmin, xmax, 1000)
        p = norm.pdf(x, mu, std )
        ax.plot(x, p, 'k', linewidth=2)
        ax.axvline(mu,**mu_kwargs)
        for num in range(1,n_sig+1):
            ax.axvline(mu+num*std,alpha=1/num,**sigma_kwargs)
            ax.axvline(mu-num*std,alpha=1/num,**sigma_kwargs)
    
    # title = "Fit results: mu = %.2f,  std = %.2f" % (mu, std)

    # print(max(p))

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    
    plt.show()
    return

def acq_parser(acq):
    if acq == '13550':
        return "Run6b_"
    elif acq == '13401':
        return "Run6_"
    else:
        return -1

def dstype_parser(id):
    if id== "dark_defects_results":
        return "DarkPixels"
    elif id== "bright_defects_results":
        return "BrightPixels"
    else: 
        return -1

def dstype_indexer(id):
    if id=="dark_defects_results":
        return "dark_pixels"
    elif id=="bright_defects_results":
        return "bright_pixels"
    else:
        return -1

def findMinMidline(arr,suggestion=3):
    return np.min(np.sort(np.abs(1-arr[int(len(arr)/2-suggestion):int(len(arr)/2+suggestion)]) * 100)[1:])

def get_serial_from_number(numbers):
    camera = LsstCam.getCamera()
    number_arr = np.array([],dtype=str)
    for number in numbers:
        number_arr = np.append(number_arr,camera.getIdMap()[number].getSerial()[:3])
    return number_arr

def getDetectorSize(numbers):
    camera = LsstCam.getCamera()
    size_arr = np.empty(len(numbers),dtype=int)
    for number in numbers:
        size_arr[number] = camera.getIdMap()[number].getBBox().getDimensions().getY() * camera.getIdMap()[number].getBBox().getDimensions().getX()
    return size_arr

def getDetectorBBox(numbers):
    # numbers is an array of detector numbers
    camera = LsstCam.getCamera()
    size_arr = []
    for number in numbers:
        size_arr.append(camera.getIdMap()[number].getBBox())
    return size_arr

def getDetectorNumber(detectorNames):
    camera = LsstCam.getCamera()
    det_nums = {det.getName():i for i, det in enumerate(camera)}
    detector_numbers = np.empty(len(detectorNames),dtype=int)
    for name in range(len(detectorNames)):
        detector_numbers[name] = det_nums[detectorNames[name]]
    return detector_numbers

def getPercent(size,six,sixb):
    '''
    size: an array with elements corresponding the number of pixels in a detector
    six: an array corresponding to the number of defects flagged in run 6
    sixb: an array corresponding to the number of defects flagged in run 6b

    returns:
    six_adjusted: an array corresponding to the percent area of defects each detector from run 6
    sixb_adjusted: an array corresponding to the percent area of defects each detector from run 6b
    diff:  an array corresponding to the differences in percent area of defects each detector from run 6 - 6b
    '''

    six_adjusted = six/size
    sixb_adjusted = sixb/size
    diff = (six-sixb)/size
    return six_adjusted*100,sixb_adjusted*100,diff*100

def getDetectors(df):
    return np.unique(np.array(df['det_name'],dtype=str))

def clearMask(display):
    for name in display.getMaskPlaneColor().keys():
        display.setMaskTransparency(transparency=0,name=name)
    return 

def formatter(display,kw=None,mask_transparency=1):
    display.setMaskPlaneColor('SATURATED',color='lime')
    display.setMaskPlaneColor('SUSPECT',color='deepskyblue')
    display.setMaskPlaneColor('INTRP',color='hotpink')
    display.setMaskPlaneColor('SAT',color='gold')
    display.setMaskPlaneColor('UNMASKEDNAN',color='darkolivegreen')
    if kw!=None:
        colors = ['yellow','cyan']
        iter = 0
        for keyword in kw:
            display.setMaskPlaneColor(keyword,color=colors[iter])
            iter+=1
    display.setMaskTransparency(mask_transparency)
    return 

def two_detector_figure(det1,det2,img_type=["flat","bias","dark"],run_number = ["13401_w_2023_24","13550_w_2023_41"],
                        base_dir = "u/lsstccs/",tight=False,save=False,path=None,figscale=4,base_dpi=400,
                        defect_datasets=["eoDarkDefects","eoBrightDefects"],
                        defect_dataset_path=["eo_dark_defects","eo_bright_defects"]):
    
    detector_nums = [det1,det2]
    
    fig = 'mtv';
    plt.close(fig);
    fig = plt.figure(fig);
    
    fig,axs = plt.subplots(len(detector_nums)*len(run_number),len(img_type),
                           figsize=[figscale*len(img_type),figscale*len(detector_nums)*len(run_number)],
                           sharex=True,sharey=True,dpi=base_dpi)
    
    fig.suptitle("Detector {var} and {var1}, {var2} and {var3}".format(var=detector_nums[0],var1=detector_nums[1],
                                                                       var2=get_serial_from_number([detector_nums[0]])[0],
                                                                       var3=get_serial_from_number([detector_nums[1]])[0]))
    collections = 'LSSTCam/raw/all,LSSTCam/calib'.split(",") # Defining collections (that are ignored by below kernel)
    butler = daf_butler.Butler('/repo/ir2',collections=collections) # Initializing a butler with the ir2 path and the above collections
    registry = butler.registry
    det_iter = 0
    for detector in detector_nums: # A new figure
        col_iter = 0
        for img in img_type: # Populating a different column
            row_iter = 0
            for run in run_number: # populating a different row
                collection = base_dir+img+"_"+run
                kwargs={
                    "datasetType": img,
                    "collections": collection,
                    "where":
                    """
                    instrument='LSSTCam' and 
                    detector = {var}
                    """.format(var=detector)
                    }
                
                    
                datasetRefs = list(registry.queryDatasets(**kwargs))
                image = butler.get(datasetRefs[0])

                
                for datasetType,defect_path in zip(defect_datasets,defect_dataset_path):
                    image.getMask().addMaskPlane(datasetType)
                    # here, add a new mask plane to add the mask.maskPixels to
                    kwargs['datasetType'] = datasetType
                    kwargs['collections'] = base_dir+defect_path+"_"+run
                    # print("Querying",datasetType,"for mask")
                    datasetRefs = list(registry.queryDatasets(**kwargs))
                    mask = butler.get(datasetRefs[0])
                    
                    mask.maskPixels(image,maskName=datasetType)
                
                # image.getMask().removeMaskPlane("BAD")
                '''
                # TESTING HERE
                for items in datasetRefs:
                    print(items)
                '''

                # image.getMask().getMaskPlaneDict

                
                # Call the the fig her
                ax = axs[row_iter+2*det_iter,col_iter]
                
                disp = afwDisplay.Display(fig)
                disp.scale('asinh', 'zscale', Q=8)
                
                formatter(disp,kw=defect_datasets)
                plt.sca(ax)
                disp.mtv(image, title="Detector "+str(detector)+", run "+run[:5]+": "+img)
    
                row_iter+=1
            col_iter +=1
        det_iter+=1
    if tight:
        fig.tight_layout()
    plt.show()
    if save:
        fig.savefig(path,dpi=base_dpi)
    fig.clear()
    plt.close(fig)
    return

def get_axes_lims(row,col,x_max,y_max,buffer):
    if row==0:
        ylim = (y_max-buffer,y_max)
    elif row==1:
        ylim = (buffer,y_max-buffer)
    elif row==2:
        ylim = (0,buffer)
    else:
        print("something wrong with row iteration!")
        return -1

    if col==0:
        xlim = (0,buffer)
    elif col==1:
        xlim = (buffer,x_max-buffer)
    elif col==2:
        xlim = (x_max-buffer,x_max)
    else:
        print("something wrong with col iteration!")
        return -1

    return xlim,ylim

def MaskBool(maskbool):
    if maskbool:
        return "WithMask"
    else:
        return "NoMask"

def makeSimpleImg(img_type,run_number,detector,repo='/repo/ir2'):
    image = getImage(img_type,run_number,detector,repo=repo)
    fig = 'mtv';
    plt.close(fig);
    fig = plt.figure(fig);
    
    fig,axs = plt.subplots()
    fig.suptitle("{var0}, detector {var1}, Run {var3}".format(var0=img_type,var1=detector,var3=run_number))

    disp = afwDisplay.Display(fig)
    disp.scale('asinh', 'zscale', Q=8)
    
    formatter(disp)
    plt.sca(axs)
    
    disp.mtv(image) 

    return fig

def edge_detector(detector=100,figscale=1,base_dpi=180,edge_buffer=20,
                  img_type="bias",run_number = "13401_w_2023_24",
                  base_dir = "u/lsstccs/",tight=False,save=False,path=None,defect_datasets=["eoDarkDefects","eoBrightDefects"],
                  defect_dataset_path=["eo_dark_defects","eo_bright_defects"],tick_sep=1,maskTransparency=1,show=True,showmask=True,run_numbah="6"):

    collections = 'LSSTCam/raw/all,LSSTCam/calib'.split(",") # Defining collections (that are ignored by below kernel)
    butler = daf_butler.Butler('/repo/ir2',collections=collections) # Initializing a butler with the ir2 path and the above collections
    registry = butler.registry

    collection = base_dir+img_type+"_"+run_number
    kwargs={
        "datasetType": img_type,
        "collections": collection,
        "where":
        """
        instrument='LSSTCam' and 
        detector = {var}
        """.format(var=detector)
        }
        
    datasetRefs = list(registry.queryDatasets(**kwargs))
    image = butler.get(datasetRefs[0])

    maxX,maxY = image.getBBox().getMax()
    if showmask:
        for datasetType,defect_path in zip(defect_datasets,defect_dataset_path):
            image.getMask().addMaskPlane(datasetType)
            kwargs['datasetType'] = datasetType
            kwargs['collections'] = base_dir+defect_path+"_"+run_number
            datasetRefs = list(registry.queryDatasets(**kwargs))
            mask = butler.get(datasetRefs[0])
            
            mask.maskPixels(image,maskName=datasetType)
    
    fig = 'mtv';
    plt.close(fig);
    fig = plt.figure(fig);
    
    fig,axs = plt.subplots(3,3,figsize=[figscale*13,figscale*10],dpi=base_dpi)
    fig.suptitle("Detector {var1}, {var2}, Run {var3}".format(var1=detector,var2=get_serial_from_number([detector])[0],var3=run_numbah))
    row_iter = 0
    for ax1 in axs:
        col_iter = 0
        for ax in ax1:
            xlim,ylim = get_axes_lims(row_iter,col_iter,maxX,maxY,edge_buffer)

            disp = afwDisplay.Display(fig)
            disp.scale('asinh', 'zscale', Q=8)

            
            formatter(disp,kw=defect_datasets,mask_transparency=maskTransparency)
            plt.sca(ax)
            
            mtv = disp.mtv(image,**{}) 
            
            disp.setMaskTransparency(maskTransparency)           
                
            if not (col_iter==1 and row_iter==1):
                ax.grid(which='major')
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            if force_tick(col_iter):
                ax.xaxis.set_major_locator(ticker.MultipleLocator(tick_sep))
            if force_tick(row_iter):
                ax.yaxis.set_major_locator(ticker.MultipleLocator(tick_sep))

            if col_iter==0:
                ax.set_ylabel("Y Pixels")
            if row_iter==2:
                ax.set_xlabel("X Pixels")
            
            ax.set_aspect('auto')  
            col_iter+=1
        row_iter+=1

    if tight:
        fig.tight_layout()
    if save:
        fig.savefig(path,dpi=base_dpi)
    if show:
        plt.show()
    fig.clear()
    plt.close(fig)
    return 

def edge_detectorPresentation(detector=100,figscale=1,base_dpi=180,edge_buffer=20,
                  img_type="bias",run_number = "13401_w_2023_24",
                  base_dir = "u/lsstccs/",tight=False,save=False,path=None,defect_datasets=["eoDarkDefects","eoBrightDefects"],
                  defect_dataset_path=["eo_dark_defects","eo_bright_defects"],tick_sep=1,maskTransparency=1,show=True,showmask=True,run_numbah="6"):

    collections = 'LSSTCam/raw/all,LSSTCam/calib'.split(",") # Defining collections (that are ignored by below kernel)
    butler = daf_butler.Butler('/repo/ir2',collections=collections) # Initializing a butler with the ir2 path and the above collections
    registry = butler.registry

    collection = base_dir+img_type+"_"+run_number
    kwargs={
        "datasetType": img_type,
        "collections": collection,
        "where":
        """
        instrument='LSSTCam' and 
        detector = {var}
        """.format(var=detector)
        }
        
    datasetRefs = list(registry.queryDatasets(**kwargs))
    image = butler.get(datasetRefs[0])

    maxX,maxY = image.getBBox().getMax()
    if showmask:
        for datasetType,defect_path in zip(defect_datasets,defect_dataset_path):
            image.getMask().addMaskPlane(datasetType)
            kwargs['datasetType'] = datasetType
            kwargs['collections'] = base_dir+defect_path+"_"+run_number
            datasetRefs = list(registry.queryDatasets(**kwargs))
            mask = butler.get(datasetRefs[0])
            
            mask.maskPixels(image,maskName=datasetType)
    
    fig = 'mtv';
    plt.close(fig);
    fig = plt.figure(fig);
    
    fig,axs = plt.subplots(3,3,figsize=[figscale*13,figscale*10],dpi=base_dpi)
    fig.suptitle("Detector {var1}, {var2}, Run {var3}".format(var1=detector,var2=get_serial_from_number([detector])[0],var3=run_numbah))
    row_iter = 0
    for ax1 in axs:
        col_iter = 0
        for ax in ax1:
            xlim,ylim = get_axes_lims(row_iter,col_iter,maxX,maxY,edge_buffer)
            return image
            ax.imshow(image.getMaskedImage())
            
            # disp = afwDisplay.Display(fig)
            # disp.scale('asinh', 'zscale', Q=8)

            
            # formatter(disp,kw=defect_datasets,mask_transparency=maskTransparency)
            # plt.sca(ax)
            
            # mtv = disp.mtv(image,**{}) 
            
            # disp.setMaskTransparency(maskTransparency)           
                
            if not (col_iter==1 and row_iter==1):
                ax.grid(which='major')
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            if force_tick(col_iter):
                ax.xaxis.set_major_locator(ticker.MultipleLocator(tick_sep))
            if force_tick(row_iter):
                ax.yaxis.set_major_locator(ticker.MultipleLocator(tick_sep))

            if col_iter==0:
                ax.set_ylabel("Y Pixels")
            if row_iter==2:
                ax.set_xlabel("X Pixels")
            
            ax.set_aspect('auto')  
            col_iter+=1
        row_iter+=1
    # fig.colorbar()
    if tight:
        fig.tight_layout()
    if save:
        fig.savefig(path,dpi=base_dpi)
    if show:
        plt.show()
    fig.clear()
    plt.close(fig)
    return 

def force_tick(iter):
    if iter==0 or iter==2:
        return True
    else:
        return False

def in_corner_or_middle(col,row):
    if (col==0 and (row==0 or row==2)):
        return True
    elif (col==1 and row==1):
        return True
    elif (col==2 and (row==0 or row==2)):
        return True
    else:
        return False


'''
THIS IS THE OLD getImage FUNCTION
'''
def getImage_old(img_type,run_number,detector,base_dir="u/lsstccs/",repo='/repo/ir2',ref_iterator=0,band='unknown'): 
    
    '''
    getImage
    
    img_type: 
    run_number: 
    detector: The detector number of the
    base_dir="u/lsstccs/": 
    repo='/repo/ir2': 
    ref_iterator=0: returns the first datasetRef in the list - useful for datasetTypes that are unique for the FP

    
    '''
    
    collections = 'LSSTCam/raw/all,LSSTCam/calib'.split(",") # Defining collections (that are ignored by below kernel)
    butler = daf_butler.Butler(repo,collections=collections) # Initializing a butler with the ir2 path and the above collections
    registry = butler.registry

    if img_type != "raw":
        
        if int(run_number[:5])>13000: # for run 5a onwards, use this scheme
            collection = base_dir+img_type+"_"+run_number
            kwargs={
                    "datasetType": img_type,
                    "collections": collection,
                    "where":
                    """
                    instrument='LSSTCam' and 
                    detector = {var} and 
                    band = {band}
                    """.format(var=detector,band=band)}
                
            refs = list(registry.queryDatasets(**kwargs))
            return butler.get(refs[ref_iterator])
        # else: # for any data that is before run 5m, use this scheme
    
    else:
        where = "exposure.science_program=run_number"
        refs = list(butler.registry.queryDatasets("raw", where=where,  bind={"run_number": run_number}, detector=detector))
        return butler.get(refs[ref_iterator])
    

def getArray(img,run,detector,repo="repo/ir2"):
    image = getImage(img,run,detector,repo=repo)
    detector_serial = get_serial_from_number([detector])[0]
    return image.image.array,detector_serial

def fixArray(arr,axis="x"):
    if axis=='y':
        return arr,axis # Take measurements along the y axis
    elif axis!='x': 
        return -1
        # throw an error
    else:
        return arr.T,axis

def computeAverageArr(arr,avg='median'):
    if avg=="mean":
        func = np.mean # compute mean
    elif avg!='median':
        print("avg function not entered correctly")
        return -1 # Throw error
    else:
        func = np.median # compute along median
    new_arr,new_arr2 = np.empty(0),np.empty(0)
    for row in arr:
        new_arr = np.append(new_arr,func(row))
        new_arr2 = np.append(new_arr2,np.std(row)/ func(row))
    return new_arr, new_arr2

def make1DFigure(arr,std_arr,axis,manufacturer,detector,tight=False,edge=False,midline=False,normed=False,left=True,figsize=[10,4],run='6',ylabel=None,ylim=None,levels=False,threshold=True,minor=False,LCA=False):
    if normed:
        y_string = 'Normalized counts [rel]'
    else:
        y_string = 'Counts [e-]'

    x = np.arange(len(arr))

    old_threshold = 0.8*np.max(arr)
    
    buffer = int(0.0125*len(arr))
    
    fig,ax = plt.subplots(figsize=figsize)

    
    if threshold:
        ax.axhline(old_threshold,ls = '--',color='red',label="80% threshold = {:.2f}".format(old_threshold))

    if axis=='x':
        xticks = np.arange(len(arr)+1,step=len(arr)/8)
        address = "serial"
        for x_val in xticks[1:-1]:
            ax.axvline(x_val,color='black',alpha=0.5,ls='-.')
    elif axis=='y':
        address = "parallel"
        xticks = np.arange(len(arr)+1,step=len(arr)/4)
        ax.axvline(np.median(xticks) -0.5,color='black',alpha=0.5,ls='-.')
        
    else:
        print("Axis defined incorrectly")
        return -1

    fig.suptitle("Averaged over {var2} axis for {var3} ({var1}) for run {var4}".format(var1=manufacturer,var2=address,var3=detector,var4=run))
    # ax.set_yscale('symlog')
    ax.set_xlabel("{var1} register position".format(var1=address))

    where="post"
    
    ekwargs = {"elinewidth":1,
           "capsize":5,
           "alpha":1}
    if edge:
        if left:
            ax.set_xticks([0,5,10,15,20,25,30,35])
            ax.set_xlim(0,35)
            ax.xaxis.set_minor_locator(MultipleLocator(1))
        else:
            ax.set_xticks(np.arange(len(x)-35,len(x)+1,step=5))
            ax.set_xlim(len(x)-35-1,len(x)-1)
            ax.xaxis.set_minor_locator(MultipleLocator(1))
            where="pre"
    elif midline:
        ax.set_xlim(np.median(x)-25,np.median(x)+25)
        ax.xaxis.set_minor_locator(MultipleLocator(1))
        where='mid'
        # Format for midline

    else:
        ax.set_xticks(xticks)
        ax.set_xlim(min(x)-buffer,max(x)+buffer)
        ax.xaxis.set_minor_locator(MultipleLocator(int(xticks[1]/4)))
        ekwargs = {"alpha":0.05}

    where='mid'
    
    ax.step(x,arr,color='black',where=where)
    
    ax.errorbar(x,arr,yerr=std_arr,ls="None",color='black',**ekwargs)

    if normed:
        ax.yaxis.set_minor_locator(MultipleLocator(0.05))

    if levels:
        ylim = (0.85,1.05) # Set the ylimits
        for line in np.arange(0.9,1,step=0.01):
            ax.axhline(line,ls = '--',color='blue') # Set a whole bunch of horizontal lines
    
    ax.set_ylabel(y_string)
    if ylim != None:
        ax.set_ylim(ylim)
    
    if minor:
        ax.grid()
        ax.yaxis.set_minor_locator(MultipleLocator(0.01))
        ax.grid(which='minor',axis="y",ls = '--',alpha=0.5)
    else:
        ax.grid()
    

    if LCA:
        if left:
            LCA_lvl = ax.get_xlim()[0] + 9
        else:
            LCA_lvl = ax.get_xlim()[1] - 9
        ax.axvline(LCA_lvl,0,2,ls='--',color='red',alpha=0.8,label="LCA-19636-A Requirement")

    if threshold:
        ax.legend(loc=8)

    if tight:
        fig.tight_layout()

    return fig
    
    

def normArr(arr,std_arr,normed,median=False):
    if median:
        arr = arr/np.median(arr)
        std_arr = std_arr/np.median(arr)
    elif normed:
        arr = arr/np.max(arr)
        std_arr = std_arr/np.max(arr)
    return arr,std_arr


def statisticsHandler(img,run,detector,axis,norm,median,repo="/repo/ir2"):
    # Get array
    arr,manu_type = getArray(img,run,detector,repo=repo)
    # Fix array
    arr,axis = fixArray(arr,axis=axis)
    # Compute averages
    avg_arr,std_arr = computeAverageArr(arr)
    # Norm array if needed
    avg_arr,std_arr = normArr(avg_arr,std_arr,norm,median=median)
    return avg_arr,std_arr,axis,manu_type

def generate1DStatistics_oneDetector(img,run,detector,axis,norm=False,tight=True,edge=False,midline=False,ylim=None,median=False,threshold=False,figsize=[10,4],left=True,runlabel="6",withLCA=False,repo='repo/ir2'):
    # Do all the averages
    avg_arr,std_arr,axis,manu_type = statisticsHandler(img,run,detector,axis,norm,median,repo=repo)
    # Make 1D figure
    return make1DFigure(avg_arr,std_arr,axis,manu_type,detector,tight=tight,normed=norm,edge=edge,midline=midline,ylim=ylim,threshold=threshold,figsize=figsize,left=left,run=runlabel,LCA=withLCA)


def getAllCals(img,run,axis='x',norm=True,edge=False,verbose=False,median=False):
    ITL_Arr,E2V_Arr = np.empty(0),np.empty(0)
    for number in range(205):
        avg_arr,std_arr,axis,manu_type = statisticsHandler(img,run,number,axis,norm,median)
        if get_serial_from_number([number])[0] == 'ITL':
            if axis=="x":
                ITL_Arr = np.append(ITL_Arr,avg_arr)
            elif axis=='y' and len(avg_arr)==4000: # Currently excluding the corner detectors - if edge==True, then ITL array will only be the edge detectors
                ITL_Arr = np.append(ITL_Arr,avg_arr)
            elif (axis=='y' and edge):
                ITL_Arr = np.append(ITL_Arr,avg_arr)
        else:
            E2V_Arr = np.append(E2V_Arr,avg_arr) 
        if verbose and number%5==0: # Every five detectors
            print("Finished detector",number)

    if axis=='x':
        print("x axis detected") 
        E2V_Arr_new = np.reshape(E2V_Arr,(-1,4096))
        ITL_Arr_new = np.reshape(ITL_Arr,(-1,4072))
    elif axis!='y':
        print("Axis defined incorrectly")
        return -1
    else:
        print("y axis detected")
        E2V_Arr_new = np.reshape(E2V_Arr,((-1,4004)))
        if edge:
            ITL_Arr_new = np.reshape(ITL_Arr,(-1,2000))
        else:
            ITL_Arr_new = np.reshape(ITL_Arr,(-1,4000))
    
    E2V_Normed,ITL_Normed = np.zeros(len(E2V_Arr_new[0])),np.zeros(len(ITL_Arr_new[0]))
    for ITLrow in ITL_Arr_new:
        ITL_Normed +=ITLrow
    for E2Vrow in E2V_Arr_new:
        E2V_Normed +=E2Vrow

    ITL_Normed_std, E2V_Normed_std = np.zeros(0),np.zeros(0)
    for row in ITL_Arr_new.T:
        ITL_Normed_std = np.append(ITL_Normed_std, np.std(row))
    for row in E2V_Arr_new.T:
        E2V_Normed_std = np.append(E2V_Normed_std, np.std(row))
    
    E2V_Normed,E2V_Normed_std = normArr(E2V_Normed,E2V_Normed_std,True,median=median)
    ITL_Normed,ITL_Normed_std = normArr(ITL_Normed,ITL_Normed_std,True,median=median)
    
    return E2V_Normed,ITL_Normed,E2V_Normed_std,ITL_Normed_std

def getSlice(side,arr,adjustment,endpoint):
    if side=='left':
        slce = (1-(arr/adjustment)[:endpoint]) * 100
    elif side=='right':
        slce = (1-(arr/adjustment)[-endpoint:]) * 100
    elif side=='mid':
        slce = np.abs(1-(arr/adjustment)[endpoint:-endpoint]) * 100
    else:
        print("side in getEdgeCutoffs defined incorrectly")
        return -1
    return slce

def getEdgeCutoffs(arr,adjustment=1,endpoint=35,side='left'): # why was this adjustment=0? that should have been a problem....
    slice = getSlice(side,arr,adjustment,endpoint)
    if type(slice)==int:
        return -1
    
    for number in np.arange(1,11):
        count = 0
        for entry in slice:
            if entry>number:
                count+=1
        print("For {var}%, there are {var2} pixels".format(var=number,var2=count))
    return

def getRecommendation(side):
    if side=="mid":
        return 4
    else:
        return 9

def getCount(slice,threshold): # for values in slice that are greater than threshold
    count=0
    for entry in slice:
        if entry>threshold:
            count+=1
    return count

def findEdgeCutoffs(arr,adjustment=1,endpoint=35,side='left',increment=1,percent_threshold = 50,sugOverride=False):
    # For use with serial axis only
    slce = getSlice(side,arr,adjustment,endpoint)
    
    if type(slce)==int:
        return slce

    if sugOverride==False:
        suggestion = getRecommendation(side) # the recommendation based on 
    else:
        suggestion = sugOverride

    
    count=0

    # normalize arr to adjacent amplifier to avoid gain mismatch
    if side=='left':
        newArr = arr/np.median(arr[0:NormRegionSize2(arr,"x")])
        return np.min((1-newArr[0:suggestion])*100)
    else:
        newArr = arr/np.median(arr[-NormRegionSize2(arr,"x")])
        return np.min((1-newArr[-suggestion:])*100)

    
    # 
    
    # while count<=suggestion:
    #     percent_threshold -= increment
    #     count=getCount(slice,percent_threshold)

def NormRegionSize2(arr,axis):
    if axis=="x":
        return int(len(arr)/8)
    elif axis=="y":
        return int(len(arr)/2)
    else:
        print("Axis in NormRegionSize defined incorrectly")
        return -1

def NormRegionSize(arr,axis):
    if axis=="serial":
        return len(arr)/8
    elif axis=="parallel":
        return len(arr)/2
    else:
        print("Axis in NormRegionSize defined incorrectly")
        return -1

def findAdjustment(arr,axis,side,buffer):
    edge_size = int(NormRegionSize(arr,axis))
    # print(buffer,edge_size)
    if side=='left':
        adjustment_region = arr[buffer:edge_size]
    elif side=="right":
        adjustment_region = arr[-edge_size:-buffer]
    else:
        print("side in findAdjustment defined incorrectly")
        return -1
    return np.median(adjustment_region)

def EdgeStatisticsManager(arr,axis,side,buffer):
    adjustment = findAdjustment(arr,axis,side,buffer)
    
    return getEdgeCutoffs(arr,adjustment=adjustment,side=side)

def heatmap(data, row_labels, col_labels, ax = None,
            cbar_kw = {}, cbarlabel = "", ylabel=None, xlabel=None, **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Parameters
    ----------
    data
        A 2D numpy array of shape (M, N).
    row_labels
        A list or array of length M with the labels for the rows.
    col_labels
        A list or array of length N with the labels for the columns.
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
        not provided, use current axes or create a new one.  Optional.
    cbar_kw
        A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
    cbarlabel
        The label for the colorbar.  Optional.
    **kwargs
        All other arguments are forwarded to `imshow`.
    """

    if not ax:
      ax = plt.gca()

    # Plot the heatmap
    im = ax.imshow(data,origin='lower', **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax = ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation = -90, va = "bottom")

    # Show all ticks and label them with the respective list entries.
    ax.set_xticks(np.arange(data.shape[1]), labels = col_labels)
    ax.set_yticks(np.arange(data.shape[0]), labels = row_labels)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top = False, bottom = True,
                   labeltop = False, labelbottom = True)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation = -30, ha = "right",
             rotation_mode = "anchor")

    # Turn spines off and create white grid.
    ax.spines[:].set_visible(False)

    ax.set_xticks(np.arange(data.shape[1] + 1) - 0.5, minor = True)
    ax.set_yticks(np.arange(data.shape[0] + 1) - 0.5, minor = True)
    ax.grid(which = "minor", color = "w", linestyle = '-', linewidth = 3)
    ax.tick_params(which =  "minor", bottom = False, left = False)

    if xlabel!=None:
        ax.set_xlabel(xlabel)
    if ylabel!=None:
        ax.set_ylabel(ylabel)

    return im, cbar

def calculateMatrix(row2,row):
    row2 = row2.T
    arr = np.empty((len(row),len(row2)))
    row_iter = 0
    for row_entry in row:
        col_iter = 0
        for col_entry in row2:
            arr[row_iter][col_iter] = row_entry + col_entry
            col_iter+=1
        row_iter+=1
    return arr

def linePlot(data,labels,lines,colors,markersize=10,figsize=[8,4],xrange=np.arange(1,11),xlims=(0.5,10.5),ylims=(-0.05,1.15),
             yrange=np.arange(-0,1.25,step=0.1), xlabel="Threshold for deviation from unity [%]",
             ylabel="Pictureframe area [%]"):
    
    fig,ax = plt.subplots(figsize=figsize)

    percent_levels = np.arange(1,11)

    iter=0
    for row in data:
        ax.plot(percent_levels,row,label=labels[iter],ls=lines[iter],color=colors[iter],marker=".",markersize=markersize)
        iter+=1
    ax.set_xticks(xrange)
    ax.set_yticks(yrange)
    ax.grid()
    ax.set_xlim(xlims)
    ax.set_ylim(ylims)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.legend()
    
    fig.show()

    return fig

def annotate_heatmap(im, data = None, valfmt="{x:.2f}",
                     textcolors = ("black", "white"),
                     threshold = None, **textkw):
    """
    A function to annotate a heatmap.

    Parameters
    ----------
    im
        The AxesImage to be labeled.
    data
        Data used to annotate.  If None, the image's data is used.  Optional.
    valfmt
        The format of the annotations inside the heatmap.  This should either
        use the string format method, e.g. "$ {x:.2f}", or be a
        `matplotlib.ticker.Formatter`.  Optional.
    textcolors
        A pair of colors.  The first is used for values below a threshold,
        the second for those above.  Optional.
    threshold
        Value in data units according to which the colors from textcolors are
        applied.  If None (the default) uses the middle of the colormap as
        separation.  Optional.
    **kwargs
        All other arguments are forwarded to each call to `text` used to create
        the text labels.
    """

    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = im.norm(data.max())/2.

    # Set default alignment to center, 
    # but allow it to be overwritten by textkw.
    kw = dict(horizontalalignment = "center",
              verticalalignment = "center")
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            kw.update(color=textcolors[int(im.norm(data[i, j]) > threshold)])
            text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
            texts.append(text)

    return texts


def twoPanelFig(observable,df,vert_scale = 1.05,ylim1=[0,1],ylim2=[0,5E5],ycol1=None,ycol2=None,
                ylabel1="Run 6-Run 6b [# pixels]",ylabel2=None,
                save=False,path=None,kwargs1 = {"color": 'blue',"ls": '-',"marker":"."},kwargs2 = {"color": 'black',"ls": '--'},
                kwargs3 = {"ls": '-',"marker":".",'alpha':0.5},colors = ['red','green']):
    # Dark, two rows, one column
    # upper panel has difference in darks
    # Lower panel has raw numbers for darks

    '''
    Inputs:
    observable: String, used for title and for querying the dataset
    Returns:
    
    '''

    if ycol1==None:
        ycol1 = observable+"Difference"
    if ycol2==None:
        ycol2=observable+"Pixels"
    if ylabel2==None:
        ylabel2="# pixels flagged as "+observable+" defect"
    
    fig,axs = plt.subplots(2,1,figsize=[10,6],sharex=True)
    
    fig.suptitle(observable)
    
    # Upper panel
    axs[0].plot(df["DetectorNumber"],df[ycol1],**kwargs1)
    axs[0].axhline(0,**kwargs2)
    axs[0].set_ylim(-vert_scale*np.max(np.abs(df[ycol1])),vert_scale*np.max(np.abs(df[ycol1])))
    axs[0].set_ylabel(ylabel1)
    axs[0].fill_between(df["DetectorNumber"], -vert_scale*np.max(np.abs(df[observable+"Difference"])),
                        vert_scale*np.max(np.abs(df[observable+"Difference"])), where=df["Manufacturer"] != "ITL",
                        color='grey', alpha=0.25, transform=axs[0].get_xaxis_transform(),label="E2V")
    axs[0].fill_between(df["DetectorNumber"], -vert_scale*np.max(np.abs(df[observable+"Difference"])),
                        vert_scale*np.max(np.abs(df[observable+"Difference"])), where=df["Manufacturer"] == "ITL",
                        color='white', alpha=0, transform=axs[0].get_xaxis_transform(),label="ITL")
    axs[0].legend()
    
    # Lower panel
    iter = 0 
    for run in ["Run6_","Run6b_"]:
        axs[1].plot(df["DetectorNumber"],df[ycol2[iter]],label=run[:-1],**kwargs3,c=colors[iter])
        iter+=1
    del iter
    axs[1].set_ylim(ylim2)
    axs[1].legend()
    axs[1].set_ylabel(ylabel2)
    axs[1].set_xlabel("Detector #")
    axs[1].set_xlim(-1,max(df['DetectorNumber'])+1)
    axs[1].fill_between(df["DetectorNumber"], -vert_scale*np.max(np.abs(df[observable+"Difference"])),
                        vert_scale*np.max(np.abs(df[observable+"Difference"])), where=df["Manufacturer"] != "ITL",
                        color='grey', alpha=0.25, transform=axs[1].get_xaxis_transform())
    
    # Final formatting
    
    for ax in axs:
        
        ax.grid(ls='--',axis='y')
    
    fig.tight_layout()
    
    plt.show()

    if save:
        fig.savefig(path+"Two_Panel.jpg",dpi=180)
    
    plt.close(fig)

    return

def getDetectors(df):
    return np.unique(np.array(df['det_name'],dtype=str))

def stackDetectors(flat_run,detectors,edge):
    all_det_avg,all_std_avg = np.array([]),np.array([])
    for detector in detectors:
        avg_arr,std_arr,__,__ = statisticsHandler('flat',flat_run,detector,'x',True,True) # x is for axis (here we want serial, so x), additional keywords are norm (bool for normalization),median (bool for normalization function)

        adjustment = findAdjustment(avg_arr,'serial',edge,20) # here, buffer = 20 and axis='serial'
        
        all_det_avg = np.append(all_det_avg,avg_arr/adjustment)
        all_std_avg = np.append(all_std_avg,std_arr/adjustment)
        
        
    all_det_avg = np.reshape(all_det_avg,(len(detectors),-1))
    all_std_avg = np.reshape(all_std_avg,(len(detectors),-1))

    all_std_avg_squared = all_std_avg*all_std_avg

    stacked_arr = np.sum(all_det_avg,axis=0)/len(detectors) # weighting each detector equally
    std_arr = np.sqrt(np.sum(all_std_avg_squared,axis=0))/np.sqrt(len(detectors)) # this is wrong...
    
    return stacked_arr,std_arr,all_det_avg,all_std_avg

def addEdgeFigData(fig,ax,stacked_data, stacked_rms, indiv_detectors,indiv_detectors_rms,detectors,edge,colors,plotIndiv=False,plotIndivErr=False):
    # Add stacked data
    if edge=='left':
        where = 'mid'
        smidge = 0
    else:
        where='mid'
        smidge=0
        
    x = np.arange(0,len(stacked_data))
    
    ax.step(x-smidge,stacked_data,color='black',marker='.',where=where,label="Stacked")
    ax.errorbar(x,stacked_data,yerr=stacked_rms,ls="None",color='black',elinewidth=1,capsize=5,alpha=1)

    if plotIndiv:
        # Add individual data
        for detector,rms,number,c in zip(indiv_detectors,indiv_detectors_rms,detectors,colors):
            ax.step(x-smidge,detector,marker='.',where=where,label=str(number),color=c,alpha=0.45)
            if plotIndivErr:
                ax.errorbar(x,detector,yerr=rms,ls="None",alpha=0.45,color=c)
    
def formatEdgeFig(fig,ax,runs,tight=False,ylims = None,ylabel=None,title=None):
    # add gridlines to all figs
    iterator=0
    for axis in ax:
        axis.grid()
    # add legend
        if iterator==0:
            axis.legend(ncols=3)
    # add labels
    # add xlim and ylim
        if ylims==None:
            axis.set_ylim(0.72,1.02)
        else:
            axis.set_ylim(ylims)
        if iterator%2==0: # on left side 
            if ylabel==None:
                axis.set_ylabel("Relative deviation [%]")
            else:
                axis.set_ylabel(ylabel)
            axis.set_xlim(-0.5,19.5)
            side='Left'
        else:
            end = 4095.5 # placeholder, will have to update if using a different raft
            axis.set_xlim(end-20,end)
            side='Right'
        if iterator<2:
            runLabel=runs[0]
        else:
            axis.set_xlabel("Serial register position")
            runLabel=runs[1]
    # Add individual titles
        if title==None:
            axis.set_title(f"{side} edge of run {runLabel}")
        else:
            axis.set_title(title[iterator])
        
        lef,rig = axis.get_xlim()
        axis.set_xticks(np.arange(int(lef),int(rig)+1,step=4))
        
        iterator+=1
    # tight layout
    if tight:
        fig.tight_layout()
    return fig

def getAmp(repo,collection,run_list,img_type,detector_num,criteria,ampName,oneRef=False):
    
    img = getImage(repo,collection,run_list,img_type,detector_num,criteria,oneRef=oneRef)
    
    detector_obj = img.getDetector()
    
    segments = [amp.getName() for amp in camera[1].getAmplifiers()]
    
    for seg in range(len(segments)):
        if segments[seg]==ampName:
            amp = detector_obj.getAmplifiers()[seg]
            break
    
    return img, detector_obj, amp

def getAmpCorners(amp):
    box = amp.getBBox()
    
    minx,maxx,miny,maxy = box.getMinX(),box.getMaxX(),box.getMinY(),box.getMaxY()

    return box,minx,maxx,miny,maxy


def getRefList(registry,kwargs):
    '''
    Function to get a list of datasetRefs, given the kwargs and a registryShim
    '''
    return list(registry.queryDatasets(**kwargs))
    
def makeButler(repo,collections):
    '''
    A function to return a butler object, given the repo and collections specified
    '''
    return daf_butler.Butler(repo,collections=collections)

def makeRegistry(butler):
    '''
    Small function to convert a Butler obj to a registry shim
    '''
    return butler.registry

def makeKwargs(run_list,img_type,detector_num,collections,instrument='LSSTCam',kwargs_where=None):
    '''
    A function to make the kwargs in a butler.registry.QueryDatasets command
    '''
    kwargs = {"datasetType":"{}".format(img_type),"collections":collections,"where":"exposure.science_program in ({}) and instrument='{}' and detector={}".format(makeRunNums(run_list),instrument,detector_num)}
    if kwargs_where!=None:
        kwargs["where"] += kwargs_where
    if img_type.__contains__("efect"):
        kwargs["where"] = " ".join(kwargs["where"].split(" ")[4:])
    return kwargs

def makeRunNums(run_list):
    '''
    A function to format the run list from [12345,67890] format to the format needed for the butler query
    '''
    return str(list(map(str, run_list)))[1:-1]

def fetchAllRefs(repo,collections,run_list,img_type,detector_num,instrument='LSSTCam',kwargs_where=None):
    '''
    A function to fetch all refs meeting the criteria of the above
    '''
    butler = makeButler(repo,collections)
    registry = makeRegistry(butler)
    kwargs = makeKwargs(run_list,img_type,detector_num,collections,instrument=instrument,kwargs_where=kwargs_where)
    return getRefList(registry,kwargs)

def refSelector(refList,criteria,butler,false_criteria={},onlyUnique=True,oneRef=False):
    '''
    A function to filter out refs within a refList, based on criteria (dict)
    '''
    
    outRefs = []
    for ref in refList:
        passCriteria = True
        exp = butler.get(ref)
        mData = exp.getInfo().getMetadata().toDict()
        for key,val in zip(criteria.keys(),criteria.values()):
            if key not in mData.keys():
                passCriteria=False
                break
            elif mData[key]!=val:
                passCriteria=False
                break
        for key,val in zip(false_criteria.keys(),false_criteria.values()):
            if mData[key]==val:
                passCriteria=False
                break
        if passCriteria:
            outRefs.append(ref)
            
        if len(outRefs)>0 and oneRef:
            return outRefs

    if onlyUnique:
        return np.unique(outRefs)
    else:
        return outRefs

def getMetadataDict(ref,butler):
    '''
    Function to get metadata of a ref as a dict
    '''
    return butler.get(ref).getInfo().getMetadata().toDict()

def getImage(repo,collection,run_list,img_type,detector_num,criteria,oneRef=False):
    '''
    A function to return an ExposureF based on criteria, repo, collection, run_list, etc....
    '''
    refs = fetchAllRefs(repo,collection,run_list,img_type,detector_num)
    filtered_refs = refSelector(refs,criteria,makeButler(repo,collection),oneRef=oneRef)
    if len(filtered_refs)!=1:
        return makeButler(repo,collection).get(filtered_refs[0]) # Probably would be good to have a better way to select, but this works for now I guess.....
    else:
        return makeButler(repo,collection).get(filtered_refs[0])


def runISR(rawimage, doLinearize=False,doDark=False,doBias=False,doFlat=False,doDefect=False,overscanFitType='MEDIAN_PER_ROW',doParallelOverscan=False,doSaturation=False,doSaturationInterpolation=False,doApplyGains=False,usePtcGains=False):
    '''
    Function to run ISR on a raw image
    '''
    isr = IsrTask()
    isr.config.doLinearize=doLinearize
    isr.config.doDark=doDark
    isr.config.doBias=doBias
    isr.config.doFlat=doFlat
    isr.config.doDefect=doDefect
    isr.config.usePtcGains=usePtcGains
    isr.config.doApplyGains=doApplyGains
    isr.config.overscan.fitType= overscanFitType
    isr.config.overscan.doParallelOverscan = doParallelOverscan
    isr.config.doSaturation=doSaturation
    isr.config.doSaturationInterpolation = doSaturationInterpolation
    # isr.config.overscan.fitType: 'AKIMA_SPLINE'
    #isr.config.overscan.fitType: 'MEAN'
    postISRImg = isr.run(rawimage)
    print(isr.assembleCcd)
    return postISRImg

def generateRefList(repo,collections,run_num,img_type,detector_num,IMG_TYPE):
    '''
    A function to generate a refList of three refs that 

    INPUTS
    -----
    repo: 
    collections: 
    run_num: 
    img_type: 
    detector_num: 
    IMG_TYPE: The type of eo image test being queried

    RETURNS
    -----
    refList: An array of three datasetRefs 
    '''
    # Get all refs
    allRefs = fetchAllRefs(repo,collections,[run_num],img_type,detector_num)
    filteredRefs = refSelector(allRefs,{"IMGTYPE":IMG_TYPE},makeButler(repo,collections))
    # Now take the filtered refs, and pull the first, last, and max expTime from the set
    refList = []
    # for ref in filteredRefs:
    #     dic = getMetadataDict(ref,makeButler(repo,collections))
    #     if dic["EXPTIME"]==15:
    #         refList.append(ref) # Add the ref with the highest expTime
    #         break# Add the first ref
    for ref in filteredRefs:
        dic = getMetadataDict(ref,makeButler(repo,collections))
        if dic["EXPTIME"]==30:
            refList.append(ref) # Add the ref with the highest expTime
            break
    for ref in filteredRefs:
        dic = getMetadataDict(ref,makeButler(repo,collections))
        if dic["EXPTIME"]==180:
            refList.append(ref) # Add the ref with the highest expTime
            break
    for ref in filteredRefs:
        dic = getMetadataDict(ref,makeButler(repo,collections))
        if dic["EXPTIME"]==360:
            refList.append(ref) # Add the ref with the highest expTime
            break
    # refList.append(filteredRefs[-1]) # add the last ref
    return refList

def getAmpFromRefs(repo,collection,ref,ampName):
    
    img = makeButler(repo,collection).get(ref)
    
    detector_obj = img.getDetector()
    
    segments = [amp.getName() for amp in camera[1].getAmplifiers()]
    
    for seg in range(len(segments)):
        if segments[seg]==ampName:
            amp = detector_obj.getAmplifiers()[seg]
            break
    
    return img, detector_obj, amp

def getFWHMFromRef(ref,col_number,repo,collections,ampName): # Need to fix this function 
    img, detector_obj, amp = getAmpFromRefs(repo,collections,ref,ampName)

    ISR_img = runISR(img).exposure

    
    
    column_arr = ISR_img.image.array.T[col_number] # Calling the specific column
    # column_arr = ISR_img[amp.getBBox()].image.array.T[amp.getBBox().maxX - col_number]

    newcolumn_arr = column_arr[amp.getBBox().minY:amp.getBBox().maxY+1]

    # print(column_arr,np.shape(column_arr))
    # print(column_arr[amp.getBBox().minY:amp.getBBox().maxY])
    
    halfMax = np.max(newcolumn_arr)/2 # Get the maximum value along the column, div by 2

    npix = len(newcolumn_arr[newcolumn_arr>halfMax]) # Counting the number of pixels that are above the half-maximum
    
    return npix,halfMax,newcolumn_arr # Return the number of pixels that exceed the half maximum, in that column

def getDiffDF(df1, df2,columns = ['BRIGHT_COLUMNS', 'BRIGHT_PIXELS','DARK_COLUMNS', 'DARK_PIXELS'][0:2],pix_threshold=10, column_threshold=1):
    '''
    Get two df's as outputs from the myu.eopipe_DictToDfz() function, and analyze differences between their columns, specified by the columns keyword
    '''
    
    # Use this slice for columns for now since the dark defects are complicated by the sensor edge effects
    
    # pix_threshold, column_threshold = 10,1 # Threshold in pixels, columns for a difference between EO runs to be flagged - data handling will need to be updated when dark defects are added
    thresholds = [column_threshold,pix_threshold]
    conditions = {} # Change to dict for easier data handling
    for col,threshold in zip(columns,thresholds):
        conditions[col] = threshold
    
    diffDF = pd.DataFrame()
    if (df1["BAY_SLOT"] == df2["BAY_SLOT"]).all() and (df1["SEGMENT"] == df2["SEGMENT"]).all():
        diffDF["BAY_SLOT"] = df2["BAY_SLOT"]
        diffDF["SEGMENT"]= df2["SEGMENT"]
        diffDF["BAY_SLOT_SEGMENT"] = df2["BAY_SLOT"] + "_" + df2["SEGMENT"]
        for column in columns:
            diffDF[column+"_DIFFERENCE"] = np.abs(df1[column] - df2[column])
        for condition_value,column in zip(conditions.values(),conditions.keys()):
            diffDF.drop(diffDF[diffDF[column+"_DIFFERENCE"]<condition_value].index,inplace=True)
    else:
        raise Exception("Imported dataframes do not match columns, data requires manipulation")

    return diffDF

def getEORunLabel(run_number):
    if run_number>=13006 and run_number<=13280:
        return "5"
    elif run_number>=13370 and run_number<=13447:
        return "6"
    elif run_number>=13479 and run_number<=13587:
        return "6b"

'''TESTING'''

def getDefectMask(repo,collections,run_list,defect_type,detector_num,instrument="LSSTCam",kwargs_where=None,iterator=0):
    refs = fetchAllRefs(repo,collections,run_list,defect_type,detector_num,instrument=instrument)
    mask = makeButler(repo,collections).get(refs[iterator])
    return mask