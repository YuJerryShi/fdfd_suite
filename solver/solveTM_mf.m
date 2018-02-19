function [Ez, Hx, Hy, omega] = solveTM_mf(L0, wvlen0, Omega, Nsb, xrange, yrange, eps_r, mod_reg, mod_phi, Jz0, Npml)
%% Input Parameters
% L0: length unit (e.g., L0 = 1e-9 for nm)
% wvlen0: input wavelength in L0
% Omega: modulation frequency in rad/s
% Nsb: number of sideband components
% xrange: [xmin xmax], range of domain in x-direction including PML
% yrange: [ymin ymax], range of domain in y-direction including PML
% eps_r: Nx-by-Ny array of relative permittivity
% Jz: Nx-by-Ny array of electric current source density
% Npml: [Nx_pml Ny_pml], number of cells in x- and y-normal PML

%% Output Parameters
% Ez, Hx, Hy: (2*Nsb+1)-by-1 cell of H- and E-field components for each
%       frequency sideband
% omega: (2*Nsb+1)-by-1 array of sideband frequencies

%% Set up the domain parameters.
eps0 = 8.854e-12 * L0;  % vacuum permittivity in farad/L0
mu0 = pi * 4e-7 * L0;  % vacuum permeability in henry/L0
c0 = 1/sqrt(eps0*mu0);  % speed of light in vacuum in L0/sec

N = size(eps_r);  % [Nx Ny]
L = [diff(xrange) diff(yrange)];  % [Lx Ly]
dL = L./N;  % [dx dy]

M = prod(N); 

omega0 = 2*pi*c0/wvlen0;  % angular frequency in rad/sec
n_sb = -Nsb : 1 : Nsb; 

omega = omega0 + Omega*n_sb; 
wvlen = 2*pi*c0./omega; 

%% Deal with the s_factor
Sxf = cell(2*Nsb+1, 1); 
Sxb = cell(2*Nsb+1, 1); 
Syf = cell(2*Nsb+1, 1); 
Syb = cell(2*Nsb+1, 1); 

for i = 1 : (2*Nsb + 1)
    [Sxf{i}, Sxb{i}, Syf{i}, Syb{i}] = S_create(L0, wvlen(i), xrange, yrange, N, Npml); 
end

%% Set up the permittivity and permeability in the domain.

vector_eps = reshape(eps0 * eps_r, M, 1); 

vec_phi = reshape(mod_phi, M, 1); 

T_eps = spdiags(vector_eps, 0, M, M); 

T_phi = spdiags(exp(1i*vec_phi), 0, M, M); 

%% Construct derivate matrices
Dyb = createDws('y', 'b', dL, N); 
Dxb = createDws('x', 'b', dL, N); 
Dxf = createDws('x', 'f', dL, N); 
Dyf = createDws('y', 'f', dL, N); 

%% Reshape Mz into a vector
jz0 = reshape(Jz0, M, 1); 

%% Construct A matrix and b vector

A_i = cell(2*Nsb+1, 1); 

for i = 1 : (2*Nsb + 1)
    A_i{i} = Sxb{i}*Dxb * mu0^-1 * Sxf{i}*Dxf + Syb{i}*Dyb * mu0^-1 * Syf{i}*Dyf + omega(i)^2*T_eps; 
end

size(omega(Nsb+1))
size(jz0)
b0 = 1i * omega(Nsb+1) * jz0;

b = zeros(M * (2*Nsb+1), 1); 

b((Nsb*M)+1 : (Nsb+1)*M, 1) = b0; 

%% Account for coupling
tic
if (Nsb > 0)
    delta_vec = reshape(eps0*mod_reg, M, 1); 
    T_delta = spdiags(delta_vec, 0, M, M);
    
    C_p = kron(spdiags([0,omega(1:end-1).^2/2]', 1, 2*Nsb+1, 2*Nsb+1), T_delta*conj(T_phi));
    C_m = kron(spdiags([omega(2:end).^2/2,0]', -1, 2*Nsb+1, 2*Nsb+1), T_delta * T_phi);
    A = blkdiag(A_i{:})+C_p+C_m;

else
    A(1:M, 1:M) = A_i{1}; 
end
time_setup = toc


%% Solve the equation.
tic
if all(b==0)
	ez = zeros(size(b));
else
	ez = A\b;
end
time_solve = toc

%%
ez_i = cell(2*Nsb+1, 1); 
hx_i = cell(2*Nsb+1, 1); 
hy_i = cell(2*Nsb+1, 1); 

Ez = cell(2*Nsb+1, 1); 
Hx = cell(2*Nsb+1, 1); 
Hy = cell(2*Nsb+1, 1); 

for i = 1 : (2*Nsb+1)
    ez_i{i} = ez((i-1)*M + 1 : i*M, 1); 
    
    hx_i{i} = -1/(1i*omega(i)) * mu0^-1 * Syf{i}*Dyf * ez_i{i}; 
    hy_i{i} = 1/(1i*omega(i)) * mu0^-1 * Sxf{i}*Dxf * ez_i{i}; 
    
    
    Ez{i} = reshape(ez_i{i}, N); 
    Hx{i} = reshape(hx_i{i}, N); 
    Hy{i} = reshape(hy_i{i}, N); 
end

end
