function [Hz, Ex, Ey, A, omega] = solveTE_modes(L0, wvlen, xrange, yrange, eps_r, Npml, n)
%% Input Parameters
% L0: length unit (e.g., L0 = 1e-9 for nm)
% wvlen: wavelength in L0
% xrange: [xmin xmax], range of domain in x-direction including PML
% yrange: [ymin ymax], range of domain in y-direction including PML
% eps_r: Nx-by-Ny array of relative permittivity
% Npml: [Nx_pml Ny_pml], number of cells in x- and y-normal PML
% n: number of eigen-modes 

%% Output Parameters
% Hz, Ex, Ey: Nx-by-Ny arrays of H- and E-field components
% A: system matrix of A x = lambda x
% omega: n-by-1 array of eigen-frequency for each mode

%% Set up the domain parameters.
eps0 = 8.854e-12 * L0;  % vacuum permittivity in farad/L0
mu0 = pi * 4e-7 * L0;  % vacuum permeability in henry/L0
c0 = 1/sqrt(eps0*mu0);  % speed of light in vacuum in L0/sec

N = size(eps_r);  % [Nx Ny]
L = [diff(xrange) diff(yrange)];  % [Lx Ly]
dL = L./N;  % [dx dy]

M = prod(N); 

omega0 = 2*pi*c0/wvlen;  % angular frequency in rad/sec

%% Deal with the s_factor
[Sxf, Sxb, Syf, Syb] = S_create(L0, wvlen, xrange, yrange, N, Npml); 

%% Set up the permittivity and permeability in the domain.
eps_x = bwdmean_w(eps0*eps_r, 'x'); 
eps_y = bwdmean_w(eps0*eps_r, 'y'); 

% Setup the Teps_x, Teps_y, and Tmu_z matrices
T_eps_x = spdiags(eps_x(:), 0, M, M); 
T_eps_y = spdiags(eps_y(:), 0, M, M); 


%% Construct derivate matrices
Dyb = Syb*createDws('y', 'b', dL, N); 
Dxb = Sxb*createDws('x', 'b', dL, N); 
Dxf = Sxf*createDws('x', 'f', dL, N); 
Dyf = Syf*createDws('y', 'f', dL, N); 

%% Construct A matrix and b vector
A = Dxf*(T_eps_x^-1)*Dxb + Dyf*(T_eps_y^-1)*Dyb; 

%% Solve the eigenvalue equation.
[hz, eig_vals] = eigs(A, n, -omega0^2*mu0); 

omega = diag(sqrt(-eig_vals/mu0)); 

%% Obtain the actual fields
Hz = cell(n, 1); 
Ex = cell(n, 1); 
Ey = cell(n, 1); 

for i = 1:n
    hz_temp = hz(1:M, i); 
    ex_temp = 1/(1i*omega(i)) * T_eps_y^-1 * Dyb * hz_temp; 
    ey_temp = 1/(1i*omega(i)) * T_eps_x^-1 * (-Dxb * hz_temp); 

    Hz{i} = reshape(hz_temp, N); 
    Ex{i} = reshape(ex_temp, N); 
    Ey{i} = reshape(ey_temp, N); 
end

end
