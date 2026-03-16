# gbt-edge-analysis

The repository includes scripts and parameters used in Teng et al. (2026, ApJ submitted) for analyzing the GBT-EDGE survey data.   

This document describes the general workflow after acquiring the reduced CO data cubes, which can be downloaded from the [EDGE team website](https://pages.astro.umd.edu/~bolatto/EDGE/).
Our data reduction pipeline is also publicly available [here](https://github.com/teuben/GBT-EDGE).

## Required Data

* CALIFA Pipe3D cubes: https://ifs.astroscu.unam.mx/CALIFA/V500/v2.3/pyPipe3D/
* Reduced CO data cubes and/or maps: https://pages.astro.umd.edu/~bolatto/EDGE/#data
* Mega-table for all the basic and derived parameters: *tables/galaxy_parameters.csv*


## Map Products

This section can be skipped if using our map products directly (see download link above).  

### 1. Pre-processing 

* Create folders to store input/output files: e.g., *data/*, *maps/*, *masks/*, *plots/*, ..., etc.
* Set up *tables/galaxy_list.csv*: define mask versions and the data sessions to be included for each galaxy
* Set up *tables/galaxy_mask.csv*: list basic parameters and masking methods to be used for each run
* Run *add_celestial.py*: add WCS info into all the CALIFA Pipe3D cubes
* Run *pipe3d_combine.py*: combine Pipe3D data for NGC0169 and NGC5929 (interacting galaxy pairs)

### 2. Masks and moment maps

* *mkmaskGBT.m*: create CALIFA-based masks for all galaxies using info from *tables/galaxy_mask.csv* 
* *mkGBTmaps.py*: define the functions needed to run *autorunGBTmaps.py* 
* *autorunGBTmaps.py*: apply masks to CO data cubes and produce resulting moment maps for all galaxies 
* *runORmasks.py*: create the combined "Hα + CO-dilated" masks and produce resulting moment maps for all galaxies 

### 3. Error estimation
* The flux rms error maps were already produced via *autorunGBTmaps.py* and *runORmaps.py*  
* *fake_source_loop.py*: produce additional error estimatation via a "fake source test" to account for baseline variations and masking uncertainties


## Scientific Analyses

This section includes analysis code for all the remaining results of the paper (Teng et al. 2026).

### 1. CALIFA-related quantities
* *sfr-from-sfh.py*: compute SFRs using star formation history based on CALIFA Pipe3D
* *fit_Zprime_gradient.py*: compute best-fit radial metallicity gradient for each galaxy using Curti+2017 calibrated 12+log(O/H) maps

### 2. CO-related inegrated quantities
* *compare-tdep_loop.py*: compute integrated CO flux, H2 mass, Hα-based SFR, and gas depletion time for each galaxy (via various α_CO prescriptions) 
* *flux-compare.py*: compare integrated fluxes between 9 overlapping galaxies from the GBT and ACA sample

### 3. Figures and visualization (scripts in *visualization/*)
* *retrieve-sdss.m*: retrieve SDSS 3-color images from NASA Atlas for all galaxies in *galaxy_list.csv*
* *tpeak-spec-smooth.py*: spectrally smooth CO data cubes and create Tpeak maps 
* *xxx-gallery.py*: generate a gallery of maps for all galaxies (sdss, tpeak, mom0, mom1, mom2)
* *sfr-aco-plots.py*: generate the figures that compare global SFRs and α_CO using various methods/prescriptions  
* *KS-plot.py*: generate the SFR-M_mol relation plots 
* *SFMS-plot.py*: generate the SFR-M_star plot, histograms for t_dep, and the ΔSFMS-t_dep relation plots

## Reference & Citation

If you use or reference any of these scripts in your work, please cite the following paper:

* Teng et al., "The EDGE–CALIFA Survey: Star Formation Efficiency and Galaxy Quenching across 62 Main Sequence, Green Valley, and Red Galaxies", 2026, submitted to *The Astrophysical Journal (ApJ)*. [[paper]](https://iopscience.iop.org/article/10.3847/1538-4357/ad10ae) 