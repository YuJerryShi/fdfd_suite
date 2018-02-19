function [Sx, Sy] = poyntingTM(Ez, Hx, Hy)
%% Input Parameters
% Ez, Hx, Hy: 2D arrays of E- and H-field components


%% Output Parameters
% Sx, Sy: 2D array of x- and y-components of Poynting vector

Ez_x = bwdmean_w(Ez, 'x');
Sx = -1/2 * real(Ez_x .* conj(Hy)); 

Ez_y = bwdmean_w(Ez, 'y'); 
Sy = 1/2 * real(Ez_y .* conj(Hx)); 

end