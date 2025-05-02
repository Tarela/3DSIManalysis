# 3DSIManalysis
3D‑SIM H3K27ac Domain Feature Analysis
# High-Resolution 3D-SIM H3K27ac Domain Analysis

This repository contains the code supporting the manuscript **“Super-resolution microscopy reveals distinct epigenetic states regulated by estrogen receptor activity.”** It provides scripts for processing raw 3D-SIM. czi files, including aligning multi-channel signals from the 3D-SIM data, extracting quantitative H3K27ac domain features, and generating statistics and publication-ready figures.

---

## Scripts
### 1. `channel_alignment.py`
- **Purpose:** Correct and align multi-channel (H3K27ac vs ER). For each Z-slides, we align the H3K27ac vs ER channels using 2D-cross-correlation and find the alignment with the best global matches. 

### 2. `3DSIM_pipeline.py`
- **Purpose:** Load `.czi` image stacks, parse metadata, detect cell nuclei via DAPI channel, and extract per-cell 3D domain features (volume, sphericity, boundary distance, channel correlations).  
- **Key Steps:**  
  1. Read CZI metadata and channels (DAPI, H3K27ac, ER)  
  2. Threshold and detect cell nuclei regions via DAPI channel  
  3. Compute 3D shape features
  4. Compute signal features 
  5. Output feature tables  

### 3. `H3K27ac_domainFeature_analysis.R`
- **Purpose:** Perform downstream statistical analysis and visualization of extracted domain features across treatment conditions in MCF7 cells (WT and YS).  
- **Key Steps:**  
  1. Load aggregated feature table (`H3K27ac_domain_features.txt`)[https://www.dropbox.com/scl/fi/wdqp7p0bzyjljw8twtr76/H3K27ac_domain_features.txt.zip?rlkey=mvrhp4xmowzcped4jg76oh2rs&dl=0]
  2. Conduct pairwise directional t-tests for hormonal treatments (E2, ED, TAM, FULV)  
  3. Generate violin plots, boxplots, and barplots for Figures and Supplementary Figures  

