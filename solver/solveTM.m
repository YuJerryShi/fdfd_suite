function [Ez, Hx, Hy, A, omega] = solveTM(L0, wvlen, xrange, yrange, eps_r, Jz, Npml)
%% Input Parameters
% L0: length unit (e.g., L0 = 1e-9 for nm)
% wvlen: wavelength in L0
% xrange: [xmin xmax], range of domain in x-direction including PML
% yrange: [ymin ymax], range of domain in y-direction including PML
% eps_r: Nx-by-Ny array of relative permittivity
% Jz: Nx-by-Ny array of electric current source density
% Npml: [Nx_pml Ny_pml], number of cells in x- and y-normal PML

%% Output Parameters
% Ez, Hx, Hy: Nx-by-Ny arrays of E- and H-field components
% A: system matrix of A x = b
% omega: angular frequency 

%% Set up the domain parameters.
eps0 = 8.854e-12 * L0;  % vacuum permittivity in farad/L0
mu0 = pi * 4e-7 * L0;  % vacuum permeability in henry/L0
c0 = 1/sqrt(eps0*mu0);  % speed of light in vacuum in L0/sec

N = size(eps_r);  % [Nx Ny]
L = [diff(xrange) diff(yrange)];  % [Lx Ly]
dL = L./N;  % [dx dy]

M = prod(N); 

omega = 2*pi*c0/wvlen;  % angular frequency in rad/sec

%% Deal with the s_factor
[Sxf, Sxb, Syf, Syb] = S_create(L0, wvlen, xrange, yrange, N, Npml); 

%% Set up the permittivity and permeability in the domain.


% Setup the Teps_x, Teps_y, and Tmu_z matrices
T_eps_z = spdiags(eps0*eps_r(:), 0, M, M); 

%% Construct derivate matrices
Dyb = Syb*createDws('y', 'b', dL, N); 
Dxb = Sxb*createDws('x', 'b', dL, N); 
Dxf = Sxf*createDws('x', 'f', dL, N); 
Dyf = Syf*createDws('y', 'f', dL, N); 

%% Reshape Mz into a vector
jz = reshape(Jz, M, 1); 

%% Construct A matrix and b vector
A = Dxf * mu0^-1 * Dxb + Dyf * mu0^-1 * Dyb + omega^2*T_eps_z; 
b = 1i * omega * jz; 


%% Solve the equation.
if all(b==0)
	ez = zeros(size(b));
else
	ez = A\b;
end
Ez = reshape(ez, N);

hx = -1/(1i*omega) * mu0^-1 * Dyb * ez; 
hy = 1/(1i*omega) * mu0^-1 * Dxb * ez; 

Hx = reshape(hx, N); 
Hy = reshape(hy, N); 

end
