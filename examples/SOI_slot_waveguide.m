clear; close all; clc; 
addpath('../class', '../flux', '../helper', '../solver', '../vis'); 

%%
fdfd = fdfd_wg_modes(); % Instantiate the solver

%% Initialize simulation properties
fdfd.L0 = 1e-6; % Set length scale to micron
fdfd.wvlen0 = 1.55; % Set wavelength
fdfd.xrange = [-1 1]; % [xmin, xmax]
fdfd.yrange = [-1 1]; % [ymin, ymax]
fdfd.N = [200 200]; % Number of cells, [Nx, Ny]

fdfd.n_modes = 4; 

%% Add permittivity blocks
fdfd.add_eps('rect', [-0.5, -0.2, 0.47, 0.4], 12); % Left Si waveguide
fdfd.add_eps('rect', [0.03, -0.2, 0.47, 0.4], 12); % Right Si waveguide
fdfd.add_eps('custom', @(x, y) y < -0.2, 2.25); % SiO2 substrate

%% Set up permittivity
fdfd.eps_setup(); 

%% Visualize the permittivity distribution
figure; 
fdfd.vis_structure(); 

%% Solve the fields
fdfd.simulate(); 

%% Visualize fields
mode = 1; 
figure; 
fdfd.visreal_wg(mode); 

%% Quiver plot of the cross-section fields
figure; 
fdfd.quiver_wg(mode); 
