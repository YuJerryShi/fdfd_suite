function [Ez, Hx, Hy, omega, A] = solveTM_modes(L0, wvlen0, xrange, yrange, eps_r, Npml, n)
%% Input Parameters
% L0: length unit (e.g., L0 = 1e-9 for nm)
% wvlen: wavelength in L0
% xrange: [xmin xmax], range of domain in x-direction including PML
% yrange: [ymin ymax], range of domain in y-direction including PML
% eps_r: Nx-by-Ny array of relative permittivity
% Npml: [Nx_pml Ny_pml], number of cells in x- and y-normal PML
% n: number of eigen-modes 

%% Output Parameters
% Ez, Hx, Hy: Nx-by-Ny arrays of H- and E-field components
% A: system matrix of A x = lambda x
% omega: n-by-1 array of eigen-frequency for each mode

%% Set up the domain parameters.
eps0 = 8.854e-12 * L0;  % vacuum permittivity in farad/L0
mu0 = pi * 4e-7 * L0;  % vacuum permeability in henry/L0
c0 = 1/sqrt(eps0*mu0);  % speed of light in vacuum in L0/sec

N = size(eps_r);  % [Nx Ny]
L = [diff(xrange) diff(yrange)];  % [Lx Ly]
dL = L./(N);  % [dx dy]

M = prod(N); 

omega0 = 2*pi*c0/wvlen0;  % angular frequency in rad/sec


%% Deal with the s_factor
[Sxf, Sxb, Syf, Syb] = S_create(L0, wvlen0, xrange, yrange, N, Npml); 


%% Set up the permittivity and permeability in the domain.
% Reshape epsilon into 1D array 
vector_eps_r = reshape(eps_r, M, 1); 

% Setup the Teps_x, Teps_y, and Tmu_z matrices
T_eps_r = spdiags(vector_eps_r, 0, M, M); 


%% Construct derivate matrices
Dyb = Syb * createDws('y', 'b', dL, N); 
Dxb = Sxb * createDws('x', 'b', dL, N); 
Dxf = Sxf * createDws('x', 'f', dL, N); 
Dyf = Syf * createDws('y', 'f', dL, N); 


%% Construct A matrix and b vector
A = T_eps_r^-1 * (Dxb * Dxf + Dyb * Dyf); 

%% Solve for modes

[ez, eig_vals] = eigs(A, n, -omega0^2*eps0*mu0); 

omega = diag(sqrt(-eig_vals/(eps0*mu0))); 

%% Obtain the actual fields
Ez = cell(n, 1); 
Hx = cell(n, 1); 
Hy = cell(n, 1); 

for i = 1:n
    ez_temp = ez(1:M, i); 
    hx_temp = -1/(1i*omega(i)) * mu0^-1 * Dyf * ez_temp; 
    hy_temp = 1/(1i*omega(i)) * mu0^-1 * Dxf * ez_temp; 
    
    Ez{i} = reshape(ez_temp, N); 
    Hx{i} = reshape(hx_temp, N); 
    Hy{i} = reshape(hy_temp, N); 
end

end

