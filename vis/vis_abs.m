function vis_abs(array2d, xrange, yrange)

%% Attach a row and column at the ends.  
% The extra row and column are not visualized, but required by pcolor().
array2d = abs(array2d);
[Nx, Ny] = size(array2d);
array2d = [array2d, array2d(:,1)];
array2d = [array2d; array2d(1,:)];

%% Create the matrices for x-locations (X), y-locations (Y), and color (C).
xs = linspace(xrange(1), xrange(2), Nx+1);
ys = linspace(yrange(1), yrange(2), Ny+1);
[X, Y] = meshgrid(xs, ys);
C = permute(array2d, [2 1]);

%% Draw with pcolor().
h = pcolor(X, Y, C);

%% Make the figure look better.
set(h, 'EdgeColor', 'none');
set(gca, 'TickDir', 'out');
axis image;

colormap('hot')
colorbar;
