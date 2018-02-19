function movie_real(array2d, xrange, yrange, Ncycle, Nfpc)

if nargin < 4  % no Ncycle
	Ncycle = 1;  % play movie for 1 cycle
end

if nargin < 5  % no Nframe
	Nfpc = 30;  % take 30 frames per cycle
end

% Multiply exp(i*phase) to array2d, and visualize the real part, where phase
% changes from 0 to 2*pi*Ncycle.
Nframe = Ncycle * Nfpc;  % total # of frames
for n = 0:Nframe
	gcf;
	vis_real(array2d * exp(1i*2*pi*n/Nfpc), xrange, yrange)
	drawnow;
end
