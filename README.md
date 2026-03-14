# gbt-edge-analysis

The repository includes scripts and parameters used in Teng et al. (2026, ApJ submitted) for analyzing the GBT-EDGE survey data.   

This document describes the general workflow after acquiring the reduced CO data, which can be downloaded from the [EDGE team website](https://pages.astro.umd.edu/~bolatto/EDGE/).
Our data reduction pipeline is also publicly available [here](https://github.com/teuben/GBT-EDGE).

## Required data products

* CALIFA Pipe3D cubes: https://ifs.astroscu.unam.mx/CALIFA/V500/v2.3/pyPipe3D/
* Reduced CO data cubes and/or maps: https://pages.astro.umd.edu/~bolatto/EDGE/#data

## Pre-processing

* Set up *galaxy_list.csv*: defines mask versions and the data sessions to be included for each galaxy
* Set up *galaxy_mask.csv*: lists basic parameters and masking methods to be used for each run
* Run *add_celestial.py*: adds WCS info into all the CALIFA Pipe3D cubes

## 1. Masks and moment maps

* *mkmaskGBT.m*: create CALIFA-based masks for all galaxies using info from *galaxy_mask.csv* 
* *mkGBTmaps.py*: define the functions needed to run *autorunGBTmaps.py* (no need to run this)
* *autorunGBTmaps.py*: apply masks to CO data cubes and produce resulting moment maps for all galaxies 
* *runORmasks.py*: create the combined "Hα + CO-dilated" masks and resulting moment maps for all galaxies 

## 2. Error estimation
* The flux rms error maps are already produced via *autorunGBTmaps.py* and *runORmaps.py*  
* *fake_source_loop.py*: produce additional error estimates via a `fake source' test to account for baseline variations and masking effects

## 3. CALIFA-related quantities
* *sfr-from-sfh.py*: compute SFRs using star formation history based on CALIFA Pipe3D
* *fit_Zprime_gradient.py*: compute best-fit radial metallicity gradient for each galaxy based on Curti et al. log(O/H) maps

## 4. CO-related inegrated quantities
* *compare-tdep_loop.py*: compute integrated CO flux, H2 mass, Hα-based SFR, and gas depletion time for each galaxy (via various α_CO prescriptions) 
* *flux-compare.py*: compare integrated flux between 9 overlapping galaxies from the GBT and ACA sample

## Figures and visualization
* *xxx-gallery.py*: generate a figure demonstrating certain type of maps (e.g. mom0, tpeak, sdss) for all the galaxies  
* *KS-plot.py*: generate the SFR-M_mol relation plots 
* *SFMS-plot.py*: generate the SFR-M* plot, histograms for t_dep, and the dSFMS-t_dep relation plots
