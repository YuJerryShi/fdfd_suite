function [Hz, Ex, Ey, omega_eigs] = solveTE_BlochX(L0, wvlen_center, xrange, yrange, eps_r, Kx, Npml, n)
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
% Hz, Ex, Ey: n-by-1 cells of the field profile of each mode
% omega_eigs: n-by-1 array of eigen-frequency for each mode

%% Set up the domain parameters.
eps0 = 8.854e-12 * L0;  % vacuum permittivity in farad/L0
mu0 = pi * 4e-7 * L0;  % vacuum permeability in henry/L0
c0 = 1/sqrt(eps0*mu0);  % speed of light in vacuum in L0/sec

N = size(eps_r);  % [Nx Ny]
L = [diff(xrange) diff(yrange)];  % [Lx Ly]
dL = L./(N);  % [dx dy]

M = prod(N); 

omega_center = 2*pi*c0/wvlen_center;  % angular frequency in rad/sec

%% Deal with the s_factor
[Sxf, Sxb, Syf, Syb] = S_create(L0, wvlen_center, xrange, yrange, N, Npml); 

%% Set up the permittivity and permeability in the domain.
eps_x = eps0 * bwdmean_w(eps_r, 'x'); 
eps_y = eps0 * bwdmean_w(eps_r, 'y'); 
eps_z = eps0 * eps_r; 
% mu_z = mu0 .* ones(N);

% Reshape epsilon into 1D array 
vector_eps_x = reshape(eps_x, M, 1); 
vector_eps_y = reshape(eps_y, M, 1); 
vector_eps_z = reshape(eps_z, M, 1); 
% vector_mu_z = reshape(mu_z, M, 1); 

% Setup the Teps_x, Teps_y, and Tmu_z matrices
T_eps_x = spdiags(vector_eps_x, 0, M, M); 
T_eps_y = spdiags(vector_eps_y, 0, M, M); 
T_eps_z = spdiags(vector_eps_z, 0, M, M); 
% T_mu_z = spdiags(vector_mu_z, 0, M, M); 

%% Construct derivate matrices
Dyb = Syb * createDws('y', 'b', dL, N); 
Dxb = Sxb * createDws('x', 'b', dL, N); 
Dxf = Sxf * createDws('x', 'f', dL, N); 
Dyf = Syf * createDws('y', 'f', dL, N); 

%% Construct the averaging matrices
Vxb = abs(createDws('x', 'b', [2 2], N)); 
Vxf = abs(createDws('x', 'f', [2 2], N)); 
Vyb = abs(createDws('y', 'b', [2 2], N)); 
Vyf = abs(createDws('y', 'f', [2 2], N)); 

%% Construct A matrix 
size(Dxf*T_eps_x^-1*Dxb)
size(Kx)

A = -mu0^-1 * (Dxf*T_eps_x^-1*Dxb - 2*1i*Kx * Vxf*T_eps_x^-1*Dxb - Kx^2*T_eps_z^-1 + Dyf*T_eps_y^-1*Dyb); 

%% Solve the eigenvalue equation.
omega_est = omega_center; 

% [vz_temp, omega_sqr] = eigs(A, n, (omega_est)^2, 'StartVector', eps_r(:)-1); 
[vz_temp, omega_sqr] = eigs(A, n, (omega_est)^2); 

omega_eigs = sqrt(diag(omega_sqr)); 

%%
Hz = cell(n, 1); 
Ex = cell(n, 1); 
Ey = cell(n, 1); 

x_vec = transpose(linspace(xrange(1), xrange(2), N(1))); 
x_space = repmat(x_vec, [N(2), 1]); 

for i = 1:n
    hz_temp = vz_temp(:, i) .* exp(-1i * Kx * x_space);
    ex_temp = 1/(1i*omega_eigs(i)) * T_eps_y^-1 * Dyb * hz_temp; 
    ey_temp = 1/(1i*omega_eigs(i)) * T_eps_x^-1 * (-Dxb * hz_temp); 
    
    Hz{i} = reshape(hz_temp, N); 
    Ex{i} = reshape(ex_temp, N); 
    Ey{i} = reshape(ey_temp, N); 
end

%%
% 
% 
% hx = -1/(1i*omega) * mu0^-1 * Dyf * ez; 
% hy = 1/(1i*omega) * mu0^-1 * Dxf * ez; 
% 
% Hx = reshape(hx, N); 
% Hy = reshape(hy, N); 

end