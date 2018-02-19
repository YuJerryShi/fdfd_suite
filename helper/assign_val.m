function val_array = assign_val(val_array, xrange, yrange, region_cond, val_fun)
%% Credits to Wonseok Shin
%% Input parameters
% val_array: 2D array of values to be updated
% xrange: [xmin xmax], range of domain in x-direction including PML
% yrange: [xmin xmax], range of domain in y-direction including PML
% region_cond: function handle with parameters (x, y) that returns true if val_array is to be update at (x, y)
% val_fun: scalar, or function handle with parameters (x, y) that returns new value of val_array at (x,y)

%% Output Parameters
N = size(val_array);

xs = linspace(xrange(1), xrange(2), N(1) + 1);  % locations of cell edges in x
ys = linspace(yrange(1), yrange(2), N(2) + 1);  % locations of cell edges in y

xs = (xs(1:end-1) + xs(2:end)) / 2;  % locations of cell centers in x
ys = (ys(1:end-1) + ys(2:end)) / 2;  % locations of cell centers in y

[X, Y] = ndgrid(xs, ys);  % X, Y are Nx-by-Ny arrays

ind = region_cond(X, Y);  % logical indices of cells satisfying condition
if isa(val_fun, 'function_handle')  % val_fun is function handle
	val = val_fun(X, Y);  % val is Nx-by-Ny array
	val_array(ind) = val(ind);  % update val_array at locations satisfying condition
else  % val_fun is scalar
	assert(isscalar(val_fun));
	val = val_fun;
	val_array(ind) = val;  % update val_array to val at locations satisfying condition
end
