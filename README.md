# Analytical N-layer Anisotropic and Eccentric MEG Forward Model

This repository provides reference MATLAB scripts accompanying the paper:

"Analytical Solution for N-layer Anisotropic and Eccentric Spherical Head Model in MEG"
The repository provides all scripts necessary to regenerate the figures reported in the manuscript, specifically Figs. 3–13.

## Purpose
The scripts are intended as a **reference implementation** of the analytical
forward model described in the paper. They allow reproduction of the main
simulation results and serve as a starting point for further development.

## Contents
- Fig3a.m  
Script demonstrating Fig.3 (a) of my manuscript for comparison of analytical four sphere model and FEM model.

- Fig3b.m  
Script demonstrating Fig.3 (b) of my manuscript for relative error between analytical four sphere model and FEM model.

- Fig4a.m  
Script demonstrating Fig.4 (a) of my manuscript for comparison anisotropic and isotropic five sphere model.

- Fig4b.m  
Script demonstrating Fig.4 (b) of my manuscript for comparison eccentric and concentric five sphere model.

- Fig5.m  
Script demonstrating Fig.5 of my manuscript for comparison between the eccentric five sphere model (with g≈0) and the concentric five sphere model.

- Fig6.m  
Script demonstrating Fig.6 of my manuscript for comparison between the analytical eccentric five sphere  model and the numerical (FEM) model.

- Fig7.m 
Script demonstrating Fig.7 of my manuscript for magnetic field for varying distances of the electric dipole the origin.

- Fig8.m 
Script demonstrating Fig.8 of my manuscript for magnetic field for different values of g.

- Fig9.m 
Script demonstrating Fig.9 of my manuscript for magnetic field graph for the anisotropic model and the isotropic model for Qz dipole.

- Fig10.m 
Script demonstrating Fig.10 of my manuscript for the relative error of magnetic field for various nmax than to nmax =75.

- Fig11.m 
Script demonstrating Fig.11 of my manuscript for localization error as a function of dipole depth for tangential and radial sources.

- Fig12.m 
Script demonstrating Fig.12 of my manuscript for effect of dipole orientation on localization error using the proposed head model.

- Fig13.m 
Script demonstrating Fig.13 of my manuscript for amplitude bias as a function of dipole depth for the proposed model and the 3 shell model.

## Sample Data
-datax.mat
The repository includes representative sample data required to reproduce the reported results.

## Requirements
- MATLAB R2019b or newer

## Notes
These scripts are provided for research and educational purposes. The authors plan to release a more complete toolbox
in future work.
