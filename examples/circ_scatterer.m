clear; close all; clc; 
addpath('../class', '../flux', '../helper', '../solver', '../vis'); 

%%
fdfd = fdfd_solve(); % Instantiate the solver

%% Initialize simulation properties
fdfd.L0 = 1e-6; % Set length scale to micron
fdfd.wvlen0 = 1.55; % Set wavelength
fdfd.xrange = [-2 2]; % [xmin, xmax]
fdfd.yrange = [-2 2]; % [ymin, ymax]
fdfd.N = [200 200]; % Number of cells, [Nx, Ny]

fdfd.Npml = [15 15]; % Number of PMLs [Npmlx, Npmly]
fdfd.pol = 'TM'; % Set the polarization

%% Add permittivity blocks
% fdfd.add_eps('rect', [-1.5, 1.1, 2, 0.1], 20); 
% fdfd.add_eps('rect', [-0.1, 0.5, 2, 0.1], 12); 
fdfd.add_eps('circ', [0, 0, 0.3], 7); % Add a circular scatterer

%% Set up the simulation
fdfd.eps_setup(); 

%% Set up the sources
% fdfd.source_setup('point', [-0.4, 0.2], 0.5); % Point source
% [beta] = fdfd.source_setup('modal', [-1.4, -0.05], {'v', 0, 30}); % Modal source

fdfd.source_setup('tfsf', [-1.2, -1, 2.5, 2], [0.5, 40]); % Total field scattered field


%% Visualize the permittivity distribution
figure; 
fdfd.vis_structure(); 

%% Perform FDFD simulation
fdfd.simulate(); 

%% Visualize a particular field
figure; 
fdfd.visreal(fdfd.Ez); 

%% Make a movie of the propagating field
figure; 
fdfd.moviereal(fdfd.Ez); 

