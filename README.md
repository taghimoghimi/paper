# Analytical N-layer Anisotropic and Eccentric MEG Forward Model

This repository provides reference MATLAB scripts accompanying the paper:

"Analytical Solution for N-layer Anisotropic and Eccentric Spherical Head Model in MEG"

## Purpose
The scripts are intended as a **reference implementation** of the analytical
forward model described in the paper. They allow reproduction of the main
simulation results and serve as a starting point for further development.

## Contents
- main_forward_example.m  
  Example script demonstrating MEG forward computation for a single dipole.

- compute_meg_forward_nlayer.m  
  Analytical forward solver using truncated spherical harmonic expansions.

- compute_coefficients_nlayer.m  
  Computation of layer-dependent coefficients based on boundary conditions
  and anisotropy parameters.

- model_parameters_example.m  
  Definition of model geometry, conductivities, eccentricity, and truncation order.

## Requirements
- MATLAB R2019b or newer

## Notes
These scripts are provided for research and educational purposes and are not
optimized for speed. The authors plan to release a more complete toolbox
in future work.
