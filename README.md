# 3DSIManalysis
3D‑SIM H3K27ac Domain Feature Analysis
# High-Resolution 3D-SIM H3K27ac Domain Analysis

This repository contains the code supporting the Nature Biotechnology (2025) study **“High-dimensional imaging using combinatorial channel multiplexing and deep learning.”** It provides a modular pipeline for processing raw 3D-SIM `.czi` volumes, aligning multi-channel data, extracting quantitative H3K27ac domain features, and generating publication-ready figures.

---

## 🔧 Scripts

### 1. `3DSIM_pipeline.py`
- **Purpose:** Load `.czi` image stacks, parse metadata, segment nuclei via DAPI channel, and extract per-cell 3D domain features (volume, sphericity, boundary distance, channel correlations).  
- **Key Steps:**  
  1. Read CZI metadata and channels (DAPI, H3K27ac, ER)  
  2. Threshold and label DAPI for cell masks  
  3. Compute 3D shape metrics via distance transforms  
  4. Output per-slide TSV feature tables and `*_metaSummary.txt` logs  

### 2. `channel_alignment.py`
- **Purpose:** Register and correct multi-channel volumes for texture and illumination artifacts, ensuring sub-pixel alignment between DAPI, H3K27ac, and ER channels.  
- **Key Steps:**  
  1. Load fixed-reference (DAPI) and moving (H3K27ac, ER) volumes  
  2. Compute spatial transforms (e.g., mutual-information based)  
  3. Apply corrections for shading and local texture  
  4. Save aligned 3D volumes for downstream extraction  

### 3. `H3K27ac_domainFeature_analysis.R`
- **Purpose:** Perform downstream statistical analysis and visualization of extracted domain features across treatment conditions in MCF7 cells (WT and YS).  
- **Key Steps:**  
  1. Load aggregated feature table (`H3K27ac_domain_features.txt`)  
  2. Conduct pairwise directional t-tests for hormonal treatments (E2, ED, TAM, FULV)  
  3. Generate violin plots, boxplots, and barplots for Figures and Supplementary Figures  
  4. Save summary P-value tables and PDF figures  

---

## ⚙️ Installation & Dependencies

```bash
# Python 3.7+ environment
pip install imageio numpy scipy opencv-python scikit-image czifile matplotlib pandas

# R 4.0+ environment
Rscript -e 'install.packages(c("gplots","ggplot2","tidyr","gridExtra"))'
