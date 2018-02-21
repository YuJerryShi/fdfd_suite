clear all; close all; clc; 

%%
fdfd = fdfd_wg_modes(); % Instantiate the solver

%% Initialize simulation properties
fdfd.L0 = 1e-6; % Set length scale to micron
fdfd.wvlen0 = 1.5; % Set wavelength
fdfd.xrange = [-0.75 0.75]; % [xmin, xmax]
fdfd.yrange = [-0.75 0.75]; % [ymin, ymax]
fdfd.N = [200 200]; % Number of cells, [Nx, Ny]

fdfd.n_modes = 1; % Number of modes
fdfd.beta_scale = 2; % Scales the initial beta_est

%% Add permittivity blocks
eps_metal = -100 - 1i*10; 
fdfd.add_eps('rect', [-0.3, -0.1, 0.27, 0.2], eps_metal); % Left metal waveguide
fdfd.add_eps('rect', [0.03, -0.1, 0.27, 0.2], eps_metal); % Right metal waveguide

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
