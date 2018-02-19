clear; close all; clc; 
addpath('../class', '../flux', '../helper', '../solver', '../vis'); 

%%
fdfd = fdfd_blochX_modes(); % Instantiate the Bloch solver

%% Initialize simulation properties
fdfd.L0 = 1e-6; % Set length scale to micron
fdfd.wvlen0 = 1.55; % Set wavelength
fdfd.xrange = [-0.25 0.25]; % [xmin, xmax]
fdfd.yrange = [-2 2]; % [ymin, ymax]
fdfd.N = [25 200]; % Number of cells, [Nx, Ny]
fdfd.Kx = 0.0; % Bloch wavevector in 1/L0

fdfd.n_modes = 20; % Number of eigenmodes 
fdfd.pol = 'TM'; % Polarization

fdfd.Npml = [0 15]; % Number of PMLs [Npmlx, Npmly]

%% Add permittivity blocks
fdfd.add_eps('rect', [-0.5, -0.2, 1, 0.4], 12); % Add a waveguide
fdfd.add_eps('circ', [0, 0.2, 0.10], 1); % Add a hole on top of the wg

%% Set up permittivity
fdfd.eps_setup(); 

%% Visualize the permittivity distribution
figure; 
fdfd.vis_structure(); 

%% Solve the fields
fdfd.simulate(); 

%% Filter the modes 
%% This is so that you can keep only the meaningful guided resonance modes
y_locs = [-0.2 0.2]; % Filter boundary
thresh = 0.15; % Energy density threshold within boundary to keep the mode

fdfd.filtering(y_locs, thresh); 

%% Get the eigen-frequencies
fdfd.omega_eigs

%% Visualize fields
mode = 2; 
figure; 
fdfd.visreal(fdfd.Ez{mode})

