function [Sx, Sy] = poyntingTE(Hz, Ex, Ey)
%% Input Parameters
% Hz, Ex, Ey: 2D arrays of H- and E-field components


%% Output Parameters
% Sx, Sy: 2D array of x- and y-components of Poynting vector

Hz_x = bwdmean_w(Hz, 'x');
Sx = 1/2 * real(Ey .* conj(Hz_x)); 

Hz_av_y = bwdmean_w(Hz, 'y'); 
Sy = -1/2 * real(Ex .* conj(Hz_av_y)); 

end