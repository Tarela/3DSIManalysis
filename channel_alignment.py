#!/usr/bin/env python3
"""
channel_alignment.py
Author: Shengen Shawn Hu
Date: May 2025

Description:
    Processes 3D-SIM .czi images to align the H3K27ac and ER channels:
      - Saves results to a TSV file per slide
"""
## library
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
    """
    Reads a .czi file and returns image arrays and metadata.
    Returns:
        meta data of input .czi 
    """
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

def shiftCor(cropdata, shiftdata,outname,slideSIZE):
    result_array = np.zeros((41, 41))  # Covering shifts from -20 to +20 for both x and y
    cropped_ER1 = cropdata[20:(slideSIZE-20), 20:(slideSIZE-20)]
    # Loop through the x and y shifts
    for dx in range(-20, 21):
        for dy in range(-20, 21):
            # Apply shifts on ER2, keeping the size the same as the cropped_ER1
            shifted_ER2 = shiftdata[20 + dx:(slideSIZE-20) + dx, 20 + dy:(slideSIZE-20) + dy]
            # Flatten matrices for Pearson correlation
            flattened_cropped_ER1 = np.log(cropped_ER1.flatten()+1)
            flattened_shifted_ER2 = np.log(shifted_ER2.flatten()+1)
            # Calculate Pearson correlation coefficient
            correlation, _ = pearsonr(flattened_cropped_ER1, flattened_shifted_ER2)
            # Store result in the array
            result_array[dx + 20, dy + 20] = correlation

    np.savetxt("%s_shift_pearsonR.txt"%outname , result_array, delimiter="\t")
    max_index = np.argmax(result_array)
    max_position = np.unravel_index(max_index, result_array.shape)
    max_x, max_y = max_position[0] - 20, max_position[1] - 20
    plt.imshow(result_array, cmap='hot', interpolation='nearest')
    plt.title('%s sqrt-signal PearsonCor'%outname)
    plt.ylabel('y-shift, min=%s, max=%s'%(round(result_array.min(),3),round(result_array.max(),3)))
    plt.xlabel('x-shift, maxPos = %s, %s'%(max_x,max_y))
    plt.xticks(np.arange(0, 41, 5), np.arange(-20, 21, 5))
    plt.yticks(np.arange(0, 41, 5), np.arange(-20, 21, 5))
    plt.show()
    return([max_x,max_y])

inname = sys.argv[1]
slideN = sys.argv[2]
outname = inname + "_z"+slideN

in_czifile = "%s.czi"%(inname)
data = czifile.CziFile(in_czifile)
czimeta = readCZImeta(data)
channel_num = czimeta[1][2]
slide_num = czimeta[1][4]
Xlim = czimeta[1][5]
Ylim = czimeta[1][6]
if Xlim != Ylim:
    print("%s: Xlim %s != Ylim %s"%(inname, Xlim,Ylim))
#for this_slide in range(slide_num):
this_slide = int(slideN)#0
H3K27ac_data = data.asarray()[0,0,0,0,this_slide,:,:,0] # red
ER_data = data.asarray()[0,0,1,0,this_slide,:,:,0] # green
DAPI_data = data.asarray()[0,0,2,0,this_slide,:,:,0] # blue

H3K27ac_data[H3K27ac_data == 65535] = 65534
ER_data[ER_data == 65535] = 65534
DAPI_data[DAPI_data == 65535] = 65534

combined_data = np.dstack((H3K27ac_data, ER_data, DAPI_data))

chan1 = H3K27ac_data#data.asarray()[0,0,0,0,this_slide,:,:,0]
chan2 = ER_data#data.asarray()[0,0,1,0,this_slide,:,:,0]
chan3 = DAPI_data#data.asarray()[0,0,2,0,this_slide,:,:,0]

fig = plt.figure(figsize=(10, 10))
fig.add_subplot(1, 1, 1)
posMaxCor_21 = shiftCor(chan2, chan1,inname+"z%s_21"%(this_slide),Xlim)
plt.savefig("%s_shift_corHeat.png"%(outname))
plt.close()
