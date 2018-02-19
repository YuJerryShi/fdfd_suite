# fdfd_suite

<img src=https://user-images.githubusercontent.com/34921690/36367280-cfbc8ae8-1506-11e8-98e1-2cc415e2aec4.png width="200" height="150" />

# Contents
This package contains a comprehensive list of two-dimensional finite-difference frequency-domain (FDFD) programs used to simulate how electromagnetic waves interact with various dielectric/metal geometries. 

It can efficiently calculate the field patterns in both TE and TM polarizations in the following simulations: 
* **Scattering:** Compute how electromagnetic waves scatter off of a dielectric/metal object
* **Resonator mode:** Compute the natural resonating frequencies of a 2D dielectric geometry
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

# More documentations to come
I will add more documentations on how to set up your own simulations from scratch soon. For now, please refer to the **examples** folder to see how the simulations are set up. 

Please let me know if you have questions! 

Also, if you would like to have additional functionalities for your simulations, I'd be happy to discuss with you to see what I can do. 

-Jerry, 2/19/2018
