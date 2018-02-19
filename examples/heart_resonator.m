clear; close all; clc; 
addpath('../class', '../flux', '../helper', '../solver', '../vis'); 

%% Instantiate the eigenmode solver
fdfd = fdfd_modes(); 

%% Initialize simulation properties
fdfd.L0 = 1e-6; % Set length scale to micron
fdfd.wvlen0 = 3; % Set wavelength
fdfd.xrange = [-1.5 1.5]; % [xmin, xmax]
fdfd.yrange = [-1.5 1.5]; % [ymin, ymax]
fdfd.N = [200 200]; % Number of cells, [Nx, Ny]

fdfd.n_modes = 2; % Number of eigenmodes 
fdfd.pol = 'TE'; % Polarization

fdfd.Npml = [15 15]; % Number of PMLs [Npmlx, Npmly]

%% Add permittivity blocks
fdfd.add_eps('circ', [-0.25, 0.1, 0.25], 12); 
fdfd.add_eps('circ', [0.25, 0.1, 0.25], 12); 
fdfd.add_eps('rect', [-0.5, -1, 1, 1.1], 1); 


vertices = [-0.5 0.1; 0.5 0.1; 0 -0.5]; % Triangle vertices
fdfd.add_eps('poly', vertices, 12); % Add a trangular resonator

%% Set up permittivity
fdfd.eps_setup(); 

%% Visualize the permittivity distribution
figure; 
fdfd.vis_structure(); 

%% Solve the fields
fdfd.simulate(); 

%% Get the eigen-frequencies
fdfd.omega_eigs
fdfd.lambda_eigs

%% Visualize fields
mode = 1; 
figure; 
fdfd.visabs(fdfd.Hz{mode})

%% Make a movie of the fields
figure; 
fdfd.moviereal(fdfd.Hz{mode})


