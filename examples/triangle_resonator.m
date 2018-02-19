clear; close all; clc; 
addpath('../class', '../flux', '../helper', '../solver', '../vis'); 

%% Instantiate the eigenmode solver
fdfd = fdfd_modes(); 

%% Initialize simulation properties
fdfd.L0 = 1e-6; % Set length scale to micron
fdfd.wvlen0 = 1.55; % Set wavelength
fdfd.xrange = [-2 2]; % [xmin, xmax]
fdfd.yrange = [-2 2]; % [ymin, ymax]
fdfd.N = [200 200]; % Number of cells, [Nx, Ny]

fdfd.n_modes = 4; % Number of eigenmodes 
fdfd.pol = 'TE'; % Polarization

fdfd.Npml = [15 15]; % Number of PMLs [Npmlx, Npmly]

%% Add permittivity blocks
vertices = [-1 -1; 0 1; 1 -1]; % Triangle vertices
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

%% Visualize fields
mode = 4; 
figure; 
fdfd.visabs(fdfd.Hz{mode})

%% Make a movie of the fields
figure; 
fdfd.moviereal(fdfd.Hz{mode})


