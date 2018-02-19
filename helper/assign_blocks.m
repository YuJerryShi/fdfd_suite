function [array_new] = assign_blocks(array_old, input_name, loc, val, xrange, yrange, N)
%ASSIGN_BLOCKS This function adds various permittivity geometries to the
%              permittivity array of the simulation domain
%% Input parameters
%   array_old: Original permittivity array
%   input_name: String with the name of the geometry
%                  choices: 'custom', 'rect', 'circ', 'poly'
%   loc: location of the geometry
%   val: value of the permittivity 
%   xrange: [xmin, xmax] of the simulation domain
%   yrange: [ymin, ymax] of the simulation domain
%   N: [Nx, Ny] number of cells in the x and y directions
% 
%% Output parameters
%   array_new: Updated permittivity array


%% Main program
input_name = lower(input_name); 
switch input_name
    %% For 'custom', the loc = boolean function handle of x and y
    case 'custom'
        array_new = assign_val(array_old, xrange, yrange, loc, val); 
                
    %% For 'rect', loc = [x_start, y_start, x_length, y_length]
    case 'rect'
        within_rect = @(x, y) x > loc(1) & x < loc(1)+loc(3) & y > loc(2) & y < loc(2)+loc(4); 
        array_new = assign_val(array_old, xrange, yrange, within_rect, val); 
        
    %% For 'circ', loc = [x_center, y_center, radius]
    case 'circ'
        r = @(x, y) sqrt((x - loc(1)).^2 + (y - loc(2)).^2); 
        within_circ = @(x, y) r(x, y) < loc(3); 
        array_new = assign_val(array_old, xrange, yrange, within_circ, val); 
    
    %% For 'poly', loc = [vertex1; vertex2; ... vertexN]; 
    case 'poly'
        array_new = array_old; 
        
        % Polygon vertex points
        xv = loc(:, 1); 
        yv = loc(:, 2); 
        
        hx = diff(xrange)/N(1); 
        hy = diff(yrange)/N(2); 
        
        xs = linspace(xrange(1), xrange(2), N(1)); 
        ys = linspace(yrange(1), yrange(2), N(2)); 
        
        % Bound the search location
        xmin_ind = floor((min(xv) - xrange(1))/hx); 
        xmax_ind = floor((max(xv) - xrange(1))/hx)+2; 
        ymin_ind = floor((min(yv) - yrange(1))/hy); 
        ymax_ind = floor((max(yv) - yrange(1))/hy)+2; 
        
        xmin_ind = max(xmin_ind, 1); 
        xmax_ind = max(xmax_ind, N(1)); 
        ymin_ind = max(ymin_ind, 1); 
        ymax_ind = max(ymax_ind, N(2)); 
        
        % Check each grid point within bound and see if it is in polygon
        for i = xmin_ind : xmax_ind
            for j = ymin_ind : ymax_ind
                in_flag = inpolygon(xs(i), ys(j), xv, yv); 
                
                if in_flag == 1
                    array_new(i, j) = val; 
                end
            end
        end
        
end

end

