function [src_profile, beta] = solve_wg_modes_1D(pol, L0, wvlen, xrange, yrange, eps_r, n)
%SOLVE_MODE Summary of this function goes here
%   Detailed explanation goes here


%% Input Parameters
% pol: polarization, either 'TE' or 'TM'); 
% L0: length unit (e.g., L0 = 1e-9 for nm)
% wvlen: wavelength in L0
% xrange: [xmin xmax], range of domain in x-direction including PML
% yrange: [ymin ymax], range of domain in y-direction including PML
% eps_r: Nx-by-Ny array of relative permittivity
% n: order of the mode needed

%% Output Parameters
% src_profile: current density distribution for the modal source
% beta: propagation constant of the n-th mode in 1/L0

%% Set up the domain parameters.
eps0 = 8.854e-12 * L0;  % vacuum permittivity in farad/L0
mu0 = pi * 4e-7 * L0;  % vacuum permeability in henry/L0
c0 = 1/sqrt(eps0*mu0);  % speed of light in vacuum in L0/sec
omega = 2*pi*c0/wvlen;  % angular frequency in rad/sec

N = size(eps_r);  % [Nx Ny]
L = [diff(xrange) diff(yrange)];  % [Lx Ly]
dL = L./(N);  % [dx dy]

M = prod(N); 


%% Set up the permittivity and permeability in the domain.

eps_x = bwdmean_w(eps0*eps_r, 'x'); 


T_eps = spdiags(eps0*eps_r(:), 0, M, M); 
T_eps_x = spdiags(eps_x(:), 0, M, M); 



%% Construct derivate matrices
Dxb = createDws('x', 'b', dL, N); 
Dxf = createDws('x', 'f', dL, N); 

%% Construct the guesses

%% Construct A matrix
switch pol
    case 'TM' 
        A = omega^2*mu0*T_eps + Dxf * Dxb; 
        
    case 'TE'
%         A = omega^2*mu0*T_eps + T_eps*Dxf*T_eps^-1*Dxb; 
        A = omega^2*mu0*T_eps + T_eps*Dxf*T_eps_x^-1*Dxb; 
end

%% Solve the equation. 
alpha = 1; 
n_diel = sqrt(max(max(real(eps_r)))); 
beta_est = abs(2*pi*n_diel / wvlen); 

[h_temp, beta_sqr] = eigs(A, n, (beta_est)^2 * alpha); 
beta = sqrt(beta_sqr(n, n)); 

%% Output the modal pattern
switch pol
    case 'TM' 
       
        ey = h_temp(:, n); 
        
        src_profile = reshape(ey, N); 
        
        
    case 'TE'
        hy = h_temp(:, n); 
        src_profile = reshape(hy, N); 
end


end

