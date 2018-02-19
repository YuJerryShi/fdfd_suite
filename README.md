# fdfd_suite

<img src=https://user-images.githubusercontent.com/34921690/36367280-cfbc8ae8-1506-11e8-98e1-2cc415e2aec4.png width="200" height="150" />

# Contents
This package contains a comprehensive list of two-dimensional finite-difference frequency-domain (FDFD) programs used to simulate how electromagnetic waves interact with various dielectric/metal geometries. 

It can efficiently calculate the field patterns in both TE and TM polarizations in the following simulations: 
* **Scattering:** Compute how electromagnetic waves scatter off of a dielectric/metal object
* **Resonator mode:** Compute the natural resonating modes and frequencies of a 2D dielectric geometry
* **Photonic crystal:** Compute the modes supported by a periodic structure
* **Modulator:** Compute the steady-state field profile of an active device whose refractive index is modulated in time
* **Waveguide mode:** Compute the modal profile of a waveguide given its 2D cross section geometry

This repository contains the following folders: 

* **class:** This folder contains the class definitions for each FDFD simulation. These classes are used to instantiate a particular simulation. 

* **examples:** This folder contains a list of examples that set up various types of FDFD simulations. The files can be used as a template for starting your own simulation. 

* **flux:** This folder contains the programs for calculating the Poynting flux given a field distribution

* **helper:** This folder contains the helper functions in the FDFD solvers, including assigning permittivity distribution, calculating various finite-difference operators, and so on. 

* **solver:** This folder contains the implementations of the FDFD algorithm

* **vis:** This folder contains the visualization programs used for plotting the permittivity distribution, field patterns, and making a movie out of the field patterns. 

# How to run
The easiest way to learn how to use this package is to take a look at the sample programs in the **example** folder. Below is a short tutorial on what you can do with the fdfd_suite. 

## Add folders to path
First, add all folders to path. 

```
clear, close all; clc; 
addpath('../class', '../flux', '../helper', '../solver', '../vis'); 
```

## Instantiate a particular solver
* Direct solve of an electromagnetic simulation with a source: fdfd = fdfd_solve(); 
* Resonator mode solver: fdfd = fdfd_modes(); 
* Photonic crystal geometry: fdfd = fdfd_blochX_modes();
* Spatiotemporal modulation of permittivity: fdfd = fdfd_mf_solve(); 
* Waveguide cross section modes: fdfd = fdfd_wg_modes(); 

## Initialize simulation properties
While most solvers have their unique properties, the properties below are common to almost all solvers and should be initialized first. 

* **fdfd.L0**: Length scale. e.g. 1e-6 = microns
* **fdfd.wvlen0**: Operating wavelength in units of L0
* **fdfd.xrange**: Simulation domain limit in the x direction in units of L0, [xmin, xmax]
* **fdfd.yrange**: Simulation domain limit in the y direction in units of L0,  [ymin, ymax]
* **fdfd.N**: Number of cells, [Nx, Ny]
* **fdfd.pol**: Polarization, either "TE" or "TM"
* **fdfd.Npml**: Number of absorbing layers at the boundaries, i.e. PMLs, [Npmlx, Npmly]

## Add permittivity blocks


# More documentations to come
I will add more documentations on how to set up your own simulations from scratch soon. For now, please refer to the **examples** folder to see how the simulations are set up. 

Please let me know if you have questions! 

Also, if you would like to have additional functionalities for your simulations, I'd be happy to discuss with you to see what I can do. 

-Jerry, 2/19/2018
