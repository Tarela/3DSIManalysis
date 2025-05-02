#!/usr/bin/env python3
"""
3DSIM_pipeline.py
Author: Shengen Shawn Hu
Date: May 2025

Description:
    Processes 3D-SIM .czi images to extract H3K27ac domain features:
      - Reads .czi metadata
      - Performs cell segmentation on DAPI channel
      - Calculates domain features (volume, sphericity, boundary distance, channel correlations)
      - Saves results to a TSV file
"""
# ========================================
# IMPORTS: core libraries for I/O, image processing,
#           numerical routines, and statistics
# ========================================
import imageio as ii
import numpy as np
import copy
import scipy.ndimage as nd
import cv2
import scipy
import scipy.stats
import matplotlib.pyplot as plt
import sys,os
import czifile
import xml.etree.ElementTree as ET
from skimage.filters import threshold_otsu, rank
from scipy.stats import spearmanr,pearsonr

def readCZImeta(input_czi):
# ----------------------------------------
# Function: readCZImeta
#   Read and parse CZI metadata XML.
# Inputs:
#   input_czi (str): path to a .czi file
# Returns:
#   XML root element for downstream metadata queries
# ----------------------------------------
    czi = input_czi
    #### fetch meta data for later output
    metadata_raw = czi.metadata()
    # Parse the XML metadata using ElementTree
    root = ET.fromstring(metadata_raw)
    # Navigate through the XML tree to find channel information
    metall = []
    for elem in root.iter('Channel'):
        channel_name = elem.get('Name')
        # print(channel_name)
        excitation_wavelength = None
        emission_wavelength = None
        detection_range = None
        color = None
        # Finding excitation wavelength
        excitation_elem = elem.find('ExcitationWavelength')
        if excitation_elem is not None:
            excitation_wavelength = str(int(float(excitation_elem.text)))
        # Finding emission wavelength
        emission_elem = elem.find('EmissionWavelength')
        if emission_elem is not None:
            emission_wavelength = str(int(float(emission_elem.text)))
        # Finding detection wavelength range
        detection_elem = elem.find('DetectionWavelength')
        if channel_name is None or excitation_wavelength is None:
            continue
        metall.append(channel_name+"::"+excitation_wavelength+"::"+emission_wavelength)
    return [metall,input_czi.asarray().shape]

def spei2D(pic):
# ----------------------------------------
# Function: spei2D
#   Compute 2D morphological metric .
# Inputs:
#   pic (2D np.ndarray): binary mask of object
# Returns:
#   2D level morphological metric ()
# ----------------------------------------
    ultima = pic + 0
    v, h = ultima.shape
    vc, hc, = nd.center_of_mass(ultima)
    o = np.ones((v, h)) + 0
    o[np.int32(np.round(vc)), np.int32(np.round(hc))] = 0
    dist = nd.distance_transform_edt(np.int32(o))
    vals = dist * ultima
    r = vals.max()
    r = np.int32(np.ceil(r))
    pixels = np.pi * (r ** 2)
    w_ei = ultima.sum()
    b_ei = pixels - w_ei
    sp = w_ei / pixels
    return sp, w_ei, b_ei, r, pixels

def spei3D(pic):
# ----------------------------------------
# Function: spei3D
#   Compute 3D morphological metric .
# Inputs:
#   pic (3D np.ndarray): binary mask of object
# Returns:
#   2D level domain shape features
# ----------------------------------------    
    ultima = pic + 0
    z,v,h = ultima.shape
    zc, vc, hc, = nd.center_of_mass(ultima)
    o = np.ones((z,v,h))+0
    o[np.int(np.round(zc)), np.int(np.round(vc)),np.int(np.round(hc))] = 0
    dist = nd.distance_transform_edt(o)
    vals = dist*ultima
    r = vals.max()
    r = np.int(np.ceil(r))
    pixels = 8*(r**3)
    w_ei = ultima.sum()
    b_ei = pixels - w_ei
    sp = w_ei/pixels
    return sp, w_ei, b_ei, r, pixels


def Shapes(pic):
# ----------------------------------------
# Function: Shapes
#   Compute shape features for a given domain .
# Inputs:
#   binary mask of object
# Returns:
#   3D level domain shape features
# ----------------------------------------
    # calculate shape metrics
    # setup
    shapes = np.zeros((1, 5))
    # obtain sp and eis
    sp, w_ei, b_ei, r, pixels = spei3D(pic)
    # obtain volume
    shapes[0][0] = sp
    shapes[0][1] = w_ei
    shapes[0][2] = b_ei
    # obtain surface area
    sa = pic.sum() - ((nd.binary_erosion(pic, iterations=1)) + 0.0).sum()
    shapes[0][3] = sa
    # obtain sphericity
    spher = ((np.pi ** (1 / 3)) * ((6 * w_ei) ** (2 / 3))) / (sa)
    shapes[0][4] = spher
    return shapes


def Intensity_percentOV(pic1, pic2, pic3, pic1_shape, pic2_shape, pic3_shape, shape):
# ----------------------------------------
# FUNCTION: Intensity_percentOV
#   Calculate average signal of channels, and correlation between H3K27ac and ER signals
# ----------------------------------------
    # pic1 = DAPI, pic2 = H3K27ac, pic3 = ER, shape = combine
    tents = np.zeros((1, 4))
    # obtain nonbackground pixels
    vals = np.where(shape.flatten() > 0.0)
    # mean intensity value
    tents[0][0] = pic1.flatten()[vals].mean()
    tents[0][1] = pic2.flatten()[vals].mean()
    tents[0][2] = pic3.flatten()[vals].mean()
    tents[0][3] = scipy.stats.pearsonr(np.sqrt(pic2.flatten()[vals]) ,np.sqrt(pic3.flatten()[vals]))[0]
    return tents

s8 = [[1, 1, 1],
      [1, 1, 1],
      [1, 1, 1]]


### read in and preprocess .czi data
inname = sys.argv[1]
IT = 2

in_czifile = "%s.czi"%(inname)
data = czifile.CziFile(in_czifile)
czimeta = readCZImeta(data)
channel_num = czimeta[1][2]
slide_num = czimeta[1][4]
Xlim = czimeta[1][5]
Ylim = czimeta[1][6]
if Xlim != Ylim:
    print("%s: Xlim %s != Ylim %s"%(inname, Xlim,Ylim))
image_data = data.asarray()
cropped_shape = (3, slide_num, Xlim - 2*20, Xlim - 2*20  )
correctMat = np.zeros(cropped_shape, dtype=image_data.dtype)

# ----------------------------------------
# step1: analyze DAPI channel and detect cell nuc regions
# ----------------------------------------
for this_slide in range(slide_num):
    # assign 3channels
    H3K27ac_data = image_data[0,0,0,0,this_slide,:,:,0] # red
    ER_data = image_data[0,0,1,0,this_slide,:,:,0] # green
    DAPI_data = image_data[0,0,2,0,this_slide,:,:,0] # blue
    H3K27ac_data[H3K27ac_data == 65535] = 65534
    ER_data[ER_data == 65535] = 65534
    DAPI_data[DAPI_data == 65535] = 65534
    # channel alignment parameter, estimated from channel_alignment.py
    shift_param_file = "%sz%s_21_shift_pearsonR.txt"%(inname,this_slide) 
    shiftCorMat = np.loadtxt(shift_param_file, delimiter="\t")
    max_index = np.argmax(shiftCorMat)
    max_position = np.unravel_index(max_index, shiftCorMat.shape)
    [dx,dy] = max_position[0] - 20, max_position[1] - 20
    # channel alignment
    slideSIZE = Xlim
    H3K27ac_usedata = H3K27ac_data[20 + dx:(slideSIZE-20) + dx, 20 + dy:(slideSIZE-20) + dy]
    ER_usedata = ER_data[20:(slideSIZE-20), 20:(slideSIZE-20)]
    DAPI_usedata = DAPI_data[20 + dx:(slideSIZE-20) + dx, 20 + dy:(slideSIZE-20) + dy]
    correctMat[0,this_slide,:,:] = H3K27ac_usedata
    correctMat[1,this_slide,:,:] = ER_usedata
    correctMat[2,this_slide,:,:] = DAPI_usedata

#### ext size by interpolating
ext_shape = (3, slide_num + (slide_num-1)*2, Xlim - 2*20, Xlim - 2*20  )
extMat = np.zeros(ext_shape, dtype=image_data.dtype)
extMat[:,0,:,:] = correctMat[:,0,:,:]
extMat[:,ext_shape[1]-1,:,:] = correctMat[:,cropped_shape[1]-1,:,:]
realSlideAnno = np.zeros(ext_shape[1:], dtype=image_data.dtype)
realSlideAnno[0,:,:] = 1
realSlideAnno[(ext_shape[1]-1),:,:] = 1

for this_slide in range(1, (ext_shape[1]-1)):
    if this_slide % 3 == 0:
        raw_slideN = int(this_slide / 3)
        extMat[:,this_slide,:,:] = correctMat[:,raw_slideN,:,:]
        realSlideAnno[this_slide,:,:] = 1
    elif this_slide % 3 == 1:
        close_raw_slideN = int((this_slide-1) / 3)
        far_raw_slideN = int((this_slide-1) / 3) + 1
        close_raw_slide = correctMat[:,close_raw_slideN,:,:]
        far_raw_slide = correctMat[:,far_raw_slideN,:,:]
        extMat[:,this_slide,:,:] = correctMat[:,close_raw_slideN,:,:] * (2/3) + correctMat[:,far_raw_slideN,:,:] * (1/3)
    elif this_slide % 3 == 2:
        far_raw_slideN = int((this_slide-2) / 3)
        close_raw_slideN = int((this_slide-2) / 3) + 1
        close_raw_slide = correctMat[:,close_raw_slideN,:,:]
        far_raw_slide = correctMat[:,far_raw_slideN,:,:]
        extMat[:,this_slide,:,:] = correctMat[:,close_raw_slideN,:,:] * (2/3) + correctMat[:,far_raw_slideN,:,:] * (1/3)

# assign the extended 3D data
H3K27ac_usedata = extMat[0,:,:,:]
ER_usedata = extMat[1,:,:,:]
DAPI_usedata = extMat[2,:,:,:]

region_cell_cb = np.zeros(DAPI_usedata.shape)
region_nonCell_cb = np.zeros(DAPI_usedata.shape)
sobel_cb = np.zeros(DAPI_usedata.shape)
boundary_pos_cb = np.empty((0, 3))

# ----------------------------------------
# step2: detect cell boundary
# ----------------------------------------
for this_slide in range(ext_shape[1]):
    DAPI_thisslide = DAPI_usedata[this_slide,:,:]
    cb_shape = (DAPI_thisslide > 500) + 0.0
    cb_shape[cb_shape>1] = 1
    b, n = nd.label(cb_shape,structure=s8)
    simp = np.hstack(b)
    locs = np.nonzero(simp)
    counts = np.bincount(simp[locs])
    vals = counts.argsort()[::-1]
    region_cell = np.zeros(b.shape)
    region_nonCell = np.ones(b.shape)
    for i in range(15):
        edata = ((b == vals[i])) + 0.0
        pixelcount = edata.sum().astype(int)
        shapedata = nd.binary_fill_holes(edata) + 0.0
        shapeEXT = nd.binary_dilation(shapedata, iterations=10)
        shapeFill = nd.binary_fill_holes(shapeEXT) + 0.0
        region_nonCell[shapeFill > 0] = 0
        shapeInfo = spei2D(shapeFill) # return sp, w_ei, b_ei, r, pixels
        ballshape = shapeInfo[0]
        radius = shapeInfo[3]
        cellPixelSize = shapeFill.sum().astype(int)
        if  (cellPixelSize >= 100000 and ballshape > 0.4):
            region_cell[shapeFill > 0] = 1
            isCell = 1
            lastCellIndex = i
            print("slide #%s, group #%s, is cell" % (this_slide,i+1))
        else:
            isCell = 0
            print("slide #%s, group #%s, is not cell" % (this_slide, i+1))
    ### cell boundary/surface
    sobelx = cv2.Sobel(region_cell, cv2.CV_64F, 1, 0, ksize=5)
    sobely = cv2.Sobel(region_cell, cv2.CV_64F, 0, 1, ksize=5)
    sobelxy = cv2.Sobel(region_cell, cv2.CV_64F, 1, 1, ksize=5)
    sobelxy_ex = np.zeros(sobelxy.shape)
    sobelxy_ex[sobelxy != 0] = 1
    boundary_pos = np.transpose(np.vstack(sobelxy_ex.nonzero()))
    # assign to combined obj
    region_cell_cb[this_slide,:,:] = region_cell
    region_nonCell_cb[this_slide,:,:] = region_nonCell
    sobel_cb[this_slide,:,:] = sobelxy_ex
    new_boundary_pos = np.hstack(( np.transpose(np.matrix(np.repeat(this_slide,boundary_pos.shape[0]))), boundary_pos))
    boundary_pos_cb = np.vstack((boundary_pos_cb, new_boundary_pos))


# ----------------------------------------
# step3: transform to binary mask (ON pixel detection)
# ----------------------------------------
# define cell nuc regions
region_cell = region_cell_cb
region_nonCell = region_nonCell_cb
sobelxy_ex3D = sobel_cb
boundary_pos = boundary_pos_cb

# separate cell-nuc signal and non-cell-nuc signal
DAPI_data_cell = DAPI_usedata * region_cell
DAPI_data_nonCell = DAPI_usedata * region_nonCell
H3K27ac_data_cell = H3K27ac_usedata * region_cell
H3K27ac_data_nonCell = H3K27ac_usedata * region_nonCell
ER_data_cell = ER_usedata * region_cell
ER_data_nonCell = ER_usedata * region_nonCell


def cell_cutoff_hist(cell_data, nonCell_data, LAB):
    # find cutoff for real signal based on the non-cell-nuc pixel signal distribution
    NCmedian = int(np.quantile(nonCell_data, 0.5))
    NC_09 = int(np.quantile(nonCell_data, 0.9))
    NC_099 = int(np.quantile(nonCell_data, 0.99))
    NC_0999 = int(np.quantile(nonCell_data, 0.999))
    NC_09999 = int(np.quantile(nonCell_data, 0.9999))
    return (NC_09, NC_099, NC_0999, NC_09999,NCmedian)


DAPI_cutoff = cell_cutoff_hist(DAPI_data_cell, DAPI_data_nonCell, "DAPI")
H3K27ac_cutoff = cell_cutoff_hist(H3K27ac_data_cell, H3K27ac_data_nonCell, "H3K27ac")
ER_cutoff = cell_cutoff_hist(ER_data_cell, ER_data_nonCell, "ER")
# use top 1/1000 (estimated from non-cell-nuc region) as cutoff
DAPI_cutoffnc0999 = DAPI_cutoff[2]
H3K27ac_cutoffnc0999 = H3K27ac_cutoff[2]
ER_cutoffnc0999 = ER_cutoff[2]
DAPI_NCmedian = DAPI_cutoff[4]
ER_NCmedian = ER_cutoff[4]
H3K27ac_NCmedian = H3K27ac_cutoff[4]

# define ON pixel
DAPI_data_cell_nc0999 = (DAPI_data_cell > DAPI_cutoffnc0999) + 0.0
H3K27ac_data_cell_nc0999 = (H3K27ac_data_cell > H3K27ac_cutoffnc0999) + 0.0
ER_data_cell_nc0999 = (ER_data_cell > ER_cutoffnc0999) + 0.0
DAPI_data_nonCell_nc0999 = (DAPI_data_nonCell > DAPI_cutoffnc0999) + 0.0
H3K27ac_data_nonCell_nc0999 = (H3K27ac_data_nonCell > H3K27ac_cutoffnc0999) + 0.0
ER_data_nonCell_nc0999 = (ER_data_nonCell > ER_cutoffnc0999) + 0.0

### get cell ID
cellIDmat, cellNum = nd.label(region_cell)
### print slides level feature/metaData
metadata_out = [inname,slide_num,Xlim, cellNum] + \
                [DAPI_cutoffnc0999,H3K27ac_cutoffnc0999,ER_cutoffnc0999] +\
                [DAPI_NCmedian,H3K27ac_NCmedian,ER_NCmedian] 

outf = open('%s_metaSummary.txt'%(inname),'w')
outf.write("\t".join(map(str,metadata_out))+"\n")
outf.close()


# ----------------------------------------
# step4: detect H3K27ac domain and extract domain-level features
# ----------------------------------------
data_shape = H3K27ac_data_cell_nc0999 
b,n = nd.label(data_shape)
simp = np.hstack(b)
locs= np.nonzero( simp )
counts=np.bincount(simp [locs] )
vals=counts.argsort()[::-1]
outf = open("%s_H3K27acDomainFeature.txt"%(inname),'w')
for i in vals:
    # exclude domains with <10 pixels
    if int(counts[i]) < 10:
        break
    print("CD number:",i)
    # domain detection
    edata = ((b==i))+0.0
    shapedata = nd.binary_fill_holes(edata)+0.0
    shapeEXT = nd.binary_dilation(shapedata, iterations=IT)
    shapeFill = nd.binary_fill_holes(shapeEXT) + 0.0
    # find domain mass center
    thispoint = nd.center_of_mass(edata)
    # distance from mass center to nearest cell nuc boundary
    minDist = np.linalg.norm(thispoint - boundary_pos, axis=1).min()
    # find associated cell nuc
    vs_cellID = np.unique(edata * cellIDmat, return_counts=False)
    related_cellID = int(max(vs_cellID))
    newll = [inname, "ACdomain%s"%(i)] + 
            # meta and pre-identified features
            [related_cellID,counts[i],shapeFill.sum() , minDist]  + 
            # shape features
            list(Shapes(shapeFill)[0]) + 
            # signal features
            list(Intensity_percentOV(pic1 = DAPI_data_cell, pic2=H3K27ac_data_cell, pic3=ER_data_cell, pic1_shape = DAPI_data_cell_nc0999,pic2_shape = H3K27ac_data_cell_nc0999, pic3_shape = ER_data_cell_nc0999, shape=shapeFill*realSlideAnno)[0])
    outf.write("\t".join(map(str,newll))+"\n")

outf.close()
print("Done")

