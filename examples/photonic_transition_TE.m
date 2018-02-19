clear; close all; clc; 
addpath('../class', '../flux', '../helper', '../solver', '../vis'); 

%% Instantiate the MF-FDFD solver
fdfd = fdfd_mf_solve(); 

%% Initialize simulation properties
fdfd.L0 = 1e-6; % Set length scale to micron
fdfd.wvlen0 = 1.55; % Set wavelength
fdfd.xrange = [0 15]; % [xmin, xmax]
fdfd.yrange = [-1 1]; % [ymin, ymax]
fdfd.N = [600 150];  % Number of cells, [Nx, Ny]
fdfd.Nsb = 1; % Number of sideband frequencies

fdfd.Npml = [15 10]; % Number of PMLs [Npmlx, Npmly] 

fdfd.pol = 'TE'; % Polarization

%% Modulation frequency
wvlen0 = 1.55;  % 1/0.6468
wvlen1 = 1.30;   % 1/0.8879
fdfd.Omega = 2*pi*fdfd.c0/fdfd.L0* (1/wvlen1 - 1/wvlen0); 


%% Linear waveguide
a = 0.4002; 
wg_lower = -a/2; 
fdfd.add_eps('rect', [-1, wg_lower, 20, a], 12.25); 

%% Modulation region
delta = 1; % Modulation amplitude
mod_x = 1.5; % Start location of the modulation
mod_L = 13.2; % range of modulation region in x
q = 3.8124; % Modulation momentum

phi_fun = @(x,y) q*x; % Modulation phase distribution as a function handle

fdfd.add_mod_amp('rect', [mod_x, -a/2, mod_L, a/2], delta); 
fdfd.add_mod_phi('rect', [mod_x, -a/2, mod_L, a/2], phi_fun); 

%% Set up permittivity and modulation profile
fdfd.eps_setup(); 
fdfd.mod_setup(); 

%% Visualize the permittivity distribution and the modulation profile
figure; 
fdfd.vis_structure(); 

figure; 
fdfd.vis_mod(); 

%% Sources
% myfdfd.source_setup('point', [-0.4, 0.2], 0.5); 
% myfdfd.source_setup('tfsf', [-1.2, -1, 2.5, 2], [0.5, 40]); 
[beta] = fdfd.source_setup('modal', [0.8, 0.00], {'v', 0, 40}); % Set up a modal source

%%
fdfd.simulate(); 

%% Inspect various sideband frequencies
figure; 
subplot(3, 1, 1); 
fdfd.visreal(fdfd.Ey{1}); 
subplot(3, 1, 2); 
fdfd.visreal(fdfd.Ey{2}); 
subplot(3, 1, 3); 
fdfd.visreal(fdfd.Ey{3}); 

%% Total fields
Ey_tot = (fdfd.Ey{1} + fdfd.Ey{2} + fdfd.Ey{3}); 
figure; 
vis_real(Ey_tot, fdfd.xrange, fdfd.yrange)
