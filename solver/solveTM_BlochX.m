function [Ez, Hx, Hy, omega_eigs] = solveTM_BlochX(L0, wvlen_min, xrange, yrange, eps_r, Kx, Npml, n)
%% Input Parameters
% L0: length unit (e.g., L0 = 1e-9 for nm)
% wvlen_center: wavelength in L0
% xrange: [xmin xmax], range of domain in x-direction including PML
% yrange: [ymin ymax], range of domain in y-direction including PML
% eps_r: Nx-by-Ny array of relative permittivity
% Kx: Bloch wavevector in 1/L0
% Npml: [Nx_pml Ny_pml], number of cells in x- and y-normal PML
% n: number of modes 

%% Output Parameters
% Ez, Hx, Hy: n-by-1 cells of the field profile of each mode
% omega_eigs: n-by-1 array of eigen-frequency for each mode

%% Set up the domain parameters.
eps0 = 8.854e-12 * L0;  % vacuum permittivity in farad/L0
mu0 = pi * 4e-7 * L0;  % vacuum permeability in henry/L0
c0 = 1/sqrt(eps0*mu0);  % speed of light in vacuum in L0/sec

N = size(eps_r);  % [Nx Ny]
L = [diff(xrange) diff(yrange)];  % [Lx Ly]
dL = L./(N);  % [dx dy]

M = prod(N); 

omega_max = 2*pi*c0/wvlen_min;  % angular frequency in rad/sec

%% Deal with the s_factor
[Sx_f, Sx_b, Sy_f, Sy_b] = S_create(L0, wvlen_min, xrange, yrange, N, Npml); 

%% Set up the permittivity and permeability in the domain.

eps_z = eps0 * eps_r; 

% Reshape epsilon into 1D array 
vector_eps_z = reshape(eps_z, M, 1); 

% Setup the Teps_x, Teps_y, and Tmu_z matrices
T_eps_z = spdiags(vector_eps_z, 0, M, M); 

%% Construct derivate matrices
Dyb = Sy_b * createDws('y', 'b', dL, N); 
Dxb = Sx_b * createDws('x', 'b', dL, N); 
Dxf = Sx_f * createDws('x', 'f', dL, N); 
Dyf = Sy_f * createDws('y', 'f', dL, N); 

%% Construct A matrix 
A = -mu0^-1 * T_eps_z^-1 * (Dxb*Dxf - 1i*Kx*Dxb - 1i*Kx*Dxf - Kx^2*speye(M) + Dyb*Dyf); 

%% Solve the eigenvalue equation.
omega_est = omega_max; 

% [uz_temp, omega_sqr] = eigs(A, n, (omega_est)^2, 'StartVector', eps_r(:)-1); 
[uz_temp, omega_sqr] = eigs(A, n, (omega_est)^2); 

omega_eigs = sqrt(diag(omega_sqr)); 

%%
Ez = cell(n, 1); 
Hx = cell(n, 1); 
Hy = cell(n, 1); 

x_vec = transpose(linspace(xrange(1), xrange(2), N(1))); 
x_space = repmat(x_vec, [N(2), 1]); 

for i = 1:n
    ez_temp = uz_temp(:, i) .* exp(-1i * Kx * x_space); 
    hx_temp = -1/(1i*omega_eigs(i)) * mu0^-1 * Dyf * ez_temp; 
    hy_temp = 1/(1i*omega_eigs(i)) * mu0^-1 * Dxf * ez_temp; 
    
    Ez{i} = reshape(ez_temp, N); 
    Hx{i} = reshape(hx_temp, N); 
    Hy{i} = reshape(hy_temp, N); 
end

%%


end