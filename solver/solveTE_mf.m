function [Hz, Ex, Ey, omega] = solveTE_mf(L0, wvlen0, Omega, Nsb, xrange, yrange, eps_r, mod_reg, mod_phi, Mz0, Npml)
%% Input Parameters
% L0: length unit (e.g., L0 = 1e-9 for nm)
% wvlen0: input wavelength in L0
% Omega: modulation frequency in rad/s
% Nsb: number of sideband components
% xrange: [xmin xmax], range of domain in x-direction including PML
% yrange: [ymin ymax], range of domain in y-direction including PML
% eps_r: Nx-by-Ny array of relative permittivity
% Mz: Nx-by-Ny array of magnetic current source density
% Npml: [Nx_pml Ny_pml], number of cells in x- and y-normal PML

%% Output Parameters
% Hz, Ex, Ey: (2*Nsb+1)-by-1 cell of H- and E-field components for each
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

eps_x = eps0 * bwdmean_w(eps_r, 'x'); 
eps_y = eps0 * bwdmean_w(eps_r, 'y'); 

T_eps_x = spdiags(eps_x(:), 0, M, M); 
T_eps_y = spdiags(eps_y(:), 0, M, M); 


%% Construct derivate matrices
Dyb = createDws('y', 'b', dL, N); 
Dxb = createDws('x', 'b', dL, N); 
Dxf = createDws('x', 'f', dL, N); 
Dyf = createDws('y', 'f', dL, N); 

%% Reshape Mz into a vector, and convert Mz to Jx, Jy
mz0 = reshape(Mz0, M, 1); 

jx0 =  Dyb * mz0; 
jy0 = -Dxb * mz0; 


%% Construct A matrix and b vector: TO DO AFTER LUNCH, 2/16/18
A_i = cell(4*Nsb+2, 1); 

for i = 1 : (2*Nsb + 1) % WRITE THIS PART!!
    A11 = -Syb{i}*Dyb * Syf{i}*Dyf - omega(i)^2*mu0*T_eps_y; 
    A12 = Syb{i}*Dyb * Sxf{i}*Dxf; 
    A21 = Sxb{i}*Dxb * Syf{i}*Dyf; 
    A22 = -Sxb{i}*Dxb * Sxf{i}*Dxf - omega(i)^2*mu0*T_eps_x; 
    
    A_i{i} = [A11  A12; 
              A21  A22]; 
end

b0 = [jx0; jy0];

b = zeros(2*M * (2*Nsb+1), 1); 

b((Nsb*2*M)+1 : (Nsb+1)*2*M, 1) = b0; 

%% Account for coupling
tic
if (Nsb > 0)
    mod_profile_p = eps0 * mod_reg .* exp(-1i*mod_phi); 
    mod_profile_m = eps0 * mod_reg .* exp(1i*mod_phi); 
    
    
    mod_p_x = bwdmean_w(mod_profile_p, 'x'); 
    mod_p_y = bwdmean_w(mod_profile_p, 'y'); 
    mod_m_x = bwdmean_w(mod_profile_m, 'x'); 
    mod_m_y = bwdmean_w(mod_profile_m, 'y'); 
    
    
    T_mod_p_x = spdiags(mod_p_x(:), 0, M, M); 
    T_mod_p_y = spdiags(mod_p_y(:), 0, M, M); 
    T_mod_m_x = spdiags(mod_m_x(:), 0, M, M); 
    T_mod_m_y = spdiags(mod_m_y(:), 0, M, M); 
    
    T_mod_p = blkdiag(T_mod_p_y, T_mod_p_x); 
    T_mod_m = blkdiag(T_mod_m_y, T_mod_m_x); 
    
    C_p = kron(spdiags([0,omega(1:end-1).^2/2]', 1, 2*Nsb+1, 2*Nsb+1), mu0*T_mod_p);
    C_m = kron(spdiags([omega(2:end).^2/2,0]', -1, 2*Nsb+1, 2*Nsb+1), mu0*T_mod_m);
    A = blkdiag(A_i{:})+C_p+C_m;

else
    A = A_i{1}; 
end


time_setup = toc

% figure; spy(A); pause; 

%% Solve the equation.
tic
if all(b==0)
	e_vec = zeros(size(b));
else
	e_vec = A\b;
end
time_solve = toc

%%
Ex = cell(2*Nsb+1, 1); 
Ey = cell(2*Nsb+1, 1); 
Hz = cell(2*Nsb+1, 1); 


for i = 1 : (2*Nsb+1)
    start_ind = (i-1)*2*M; 
    ex_temp = e_vec(start_ind+1 : start_ind+M, 1); 
    ey_temp = e_vec(start_ind+M+1 : start_ind+2*M, 1); 
    
    hz_temp = -1/(1i*omega(i)*mu0) * (Sxf{i}*Dxf * ey_temp - Syf{i}*Dyf * ex_temp); 
    
    Ex{i} = reshape(ex_temp, N); 
    Ey{i} = reshape(ey_temp, N); 
    Hz{i} = reshape(hz_temp, N); 
end

end
