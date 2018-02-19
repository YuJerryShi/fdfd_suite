function [Sxf, Sxb, Syf, Syb] = S_create(L0, wvlen, xrange, yrange, N, Npml)
%S_CREATE: creates the stretched-coordinate PML operators in FDFD
%   Detailed explanation goes here


%% Set up the domain parameters.
eps0 = 8.854e-12 * L0;  % vacuum permittivity in farad/L0
mu0 = pi * 4e-7 * L0;  % vacuum permeability in henry/L0
c0 = 1/sqrt(eps0*mu0);  % speed of light in vacuum in L0/sec

M = prod(N); 

omega = 2*pi*c0/wvlen;  % angular frequency in rad/sec

%% Deal with the s_factor
% Create the sfactor in each direction and for 'f' and 'b'
s_vector_x_f = create_sfactor(xrange, 'f', omega, eps0, mu0, N(1), Npml(1)); 
s_vector_x_b = create_sfactor(xrange, 'b', omega, eps0, mu0, N(1), Npml(1)); 
s_vector_y_f = create_sfactor(yrange, 'f', omega, eps0, mu0, N(2), Npml(2)); 
s_vector_y_b = create_sfactor(yrange, 'b', omega, eps0, mu0, N(2), Npml(2)); 


% Fill the 2D space with layers of appropriate s-factors
Sx_f_2D = zeros(N); 
Sx_b_2D = zeros(N); 
Sy_f_2D = zeros(N); 
Sy_b_2D = zeros(N); 

for j = 1:N(2)
    Sx_f_2D(:, j) = s_vector_x_f .^-1;  
    Sx_b_2D(:, j) = s_vector_x_b .^-1; 
end

for i = 1:N(1)
    Sy_f_2D(i, :) = s_vector_y_f .^-1; 
    Sy_b_2D(i, :) = s_vector_y_b .^-1; 
end

% surf(abs(Sy_f_2D)); pause

% Reshape the 2D s-factors into a 1D s-array
Sx_f_vec = reshape(Sx_f_2D, M, 1); 
Sx_b_vec = reshape(Sx_b_2D, M, 1); 
Sy_f_vec = reshape(Sy_f_2D, M, 1); 
Sy_b_vec = reshape(Sy_b_2D, M, 1); 

% Construct the 1D total s-array into a diagonal matrix
Sxf = spdiags(Sx_f_vec, 0, M, M); 
Sxb = spdiags(Sx_b_vec, 0, M, M); 
Syf = spdiags(Sy_f_vec, 0, M, M); 
Syb = spdiags(Sy_b_vec, 0, M, M); 


end

