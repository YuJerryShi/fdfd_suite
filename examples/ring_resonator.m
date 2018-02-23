clear; close all; clc; 
addpath('../class', '../flux', '../helper', '../solver', '../vis'); 

%%
fdfd = fdfd_solve(); % Instantiate the solver

%% Initialize simulation properties
fdfd.L0 = 1e-6; % Set length scale to micron
fdfd.wvlen0 = 1.508465; % Set wavelength
fdfd.xrange = [0 6]; % [xmin, xmax]
fdfd.yrange = [0 6]; % [ymin, ymax]
fdfd.N = [300 300]; % Number of cells, [Nx, Ny]

fdfd.Npml = [15 15]; % Number of PMLs [Npmlx, Npmly]
fdfd.pol = 'TM'; % Set the polarization

%% Add permittivity blocks
r1 = 1.5; % Inner radius
r2 = 1.8; % Outer radius

wg_c = 0.81; % Waveguide center
wg_d = 0.201; % Waveguide width

fdfd.add_eps('circ', [3, 3, r2], 12); 
fdfd.add_eps('circ', [3, 3, r1], 1); 
fdfd.add_eps('rect', [fdfd.xrange(1), wg_c-wg_d/2, diff(fdfd.xrange), wg_d], 12); 

%% Set up the simulation
fdfd.eps_setup(); 

%% Set up the sources
[beta] = fdfd.source_setup('modal', [1.01, wg_c], {'v', 0, 30}); % Modal source

%% Visualize the permittivity distribution
figure; 
fdfd.vis_structure(); 

%% Perform FDFD simulation
fdfd.simulate(); 

%% Visualize a particular field
figure; 
fdfd.visabs(fdfd.Ez); 

%% Make a movie of the propagating field
figure; 
fdfd.moviereal(fdfd.Ez); 

