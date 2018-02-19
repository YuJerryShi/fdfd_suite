function [Hx, Hy, Hz, Ex, Ey, Ez, beta] = solve_wg_modes(L0, wvlen, xrange, yrange, eps_r, n)
%SOLVE_WG_MODE: Solves the waveguide mode given its x-y cross section

%% Input Parameters
% L0: length unit (e.g., L0 = 1e-6 for microns)
% wvlen: wavelength in L0
% xrange: [xmin xmax], range of domain in x-direction including PML
% yrange: [ymin ymax], range of domain in y-direction including PML
% eps_r: Nx-by-Ny array of relative permittivity
% n: number of modes

%% Output Parameters
% Hx, Hy, Hz, Ex, Ey, Ez: n-by-1 cells containing of H- and E-field
%       components of each mode
% beta: propagation constant of each mode in 1/L0

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
eps_x = bwdmean_w(eps0 * eps_r, 'y');  % average eps for eps_x
eps_y = bwdmean_w(eps0 * eps_r, 'x');  % average eps for eps_y
eps_z = bwdmean_w(eps0 * eps_r, 'x'); 
eps_z = bwdmean_w(eps_z, 'y'); % average eps for eps_z

% Reshape epsilon into 1D array 
vector_eps_x = reshape(eps_x, M, 1); 
vector_eps_y = reshape(eps_y, M, 1); 
vector_eps_z = reshape(eps_z, M, 1); 

% Setup the Teps_x, Teps_y, and Tmu_z matrices
T_eps_x = spdiags(vector_eps_x, 0, M, M); 
T_eps_y = spdiags(vector_eps_y, 0, M, M); 
T_eps_z = spdiags(vector_eps_z, 0, M, M); 

T_eps_z_inv = spdiags(vector_eps_z.^(-1), 0, M, M); 

%% Construct derivate matrices
Dyb = createDws('y', 'b', dL, N); 
Dxb = createDws('x', 'b', dL, N); 
Dxf = createDws('x', 'f', dL, N); 
Dyf = createDws('y', 'f', dL, N); 

%% Construct A matrix
A = omega^2 * mu0 * blkdiag(T_eps_y, T_eps_x) + ... 
    + blkdiag(T_eps_y, T_eps_x)* [-Dyf; Dxf] * T_eps_z_inv * [-Dyb Dxb] + ...
    + [Dxb; Dyb]*[Dxf Dyf]; 


%% Solve the equation. 
alpha = 1; 
n_diel = sqrt(max(max(real(eps_r)))); 
beta_est = abs(2*pi*n_diel / wvlen); 

[h_temp, beta_sqr] = eigs(A, n, (beta_est*2)^2 * alpha); 


%% Rearrange all the fields back into an array
Hx = cell(n, 1); 
Hy = cell(n, 1); 
Hz = cell(n, 1); 

Ex = cell(n, 1); 
Ey = cell(n, 1); 
Ez = cell(n, 1); 

beta = zeros(n, 1); 

for i = 1:n
    beta(i) = sqrt(beta_sqr(i, i)); 
    gamma = 1i*beta(i); 
    
    hx = h_temp(1:M, i); 
    hy = h_temp(M+1:2*M, i); 
    hz = 1./(gamma) * (Dxf * hx + Dyf * hy); 
    
    ex = 1./(1i*omega) * T_eps_x^-1 * (Dyb * hz + gamma * hy); 
    ey = 1./(1i*omega) * T_eps_y^-1 * (-gamma * hx - Dxb * hz); 
    ez = 1./(1i*omega) * T_eps_z^-1 * (Dxb * hy - Dyb * hx); 
    
    %% Reshape all vectors
    Hx{i} = reshape(hx, N); 
    Hy{i} = reshape(hy, N);
    Hz{i} = reshape(hz, N); 
    Ex{i} = reshape(ex, N); 
    Ey{i} = reshape(ey, N); 
    Ez{i} = reshape(ez, N); 
    
end


end

