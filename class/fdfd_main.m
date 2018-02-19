classdef fdfd_main < handle
    %fdfd_main class is a general fdfd class that contains parameters that are
    %used in more specific FDFD simulations. All other FDFD simulations
    %inherit their properties from the fdfd class
    
    
    properties (Constant)
        eps0 = 8.854e-12;  % vacuum permittivity in farad/m
        mu0 = pi * 4e-7;  % vacuum permeability in henry/m
        c0 = 1/sqrt(8.854e-12 * pi * 4e-7);  % speed of light in vacuum in m/sec
    end
    
    properties
        %% Simlulation properties
        pol % Polarization: 'TE' or 'TM'
        L0 % Length scale. 1e-6 corresponds to microns
        wvlen0 % Input wavelength in L0
        omega0 % Input frequency
        xrange % Simulation domain boundaries in L0: [xmin, xmax] 
        yrange % Simulation domain boundaries in L0: [ymin, ymax]
        N % Number of cells in x and y: [Nx, Ny]
        dL % Resolution in x and y: [dx, dy]
        Npml % Number of PMLs in x and y: [Npmlx, Npmly]
        blocks % Storage of various permittivity blocks that are added 
        eps_r % Relative permittivity distribution
        src % Source distribution of the simulation
        
        %% Output properties
        % Electric fields
        Ex 
        Ey
        Ez
        
        % Magnetic fields
        Hx
        Hy
        Hz
    end
    
    
    methods
        %% Add dielectric blocks to the simulation
        function obj = add_eps(obj, block_name, block_loc, block_val)
            N_blocks = length(obj.blocks); 

            obj.blocks{N_blocks+1} = {block_name, block_loc, block_val}; 
        end
        
        %% Visualize the permittivity distribution
        function [] = vis_structure(obj)
            vis_abs(obj.eps_r, obj.xrange, obj.yrange); 
        end
        
        %% Visualization of the real part of the field
        function [] = visreal(obj, field)
            vis_real(field, obj.xrange, obj.yrange); 
        end
        
        %% Visualization of the amplitude of the field
        function [] = visabs(obj, field)
            vis_abs(field, obj.xrange, obj.yrange); 
        end
        
        %% Make a movie out of the field
        function [] = moviereal(obj, field, Ncycle, Nfpc)
            if nargin < 4  % no Ncycle
                Ncycle = 1;  % play movie for 1 cycle
            end
            
            if nargin < 5  % no Nframe
                Nfpc = 30;  % take 30 frames per cycle
            end
            movie_real(field, obj.xrange, obj.yrange, Ncycle, Nfpc); 
        end
        
        %% Compute the Poynting flux
        function [Sx, Sy] = poynting(obj)
            switch obj.pol
                case 'TE'
                    [Sx, Sy] = poyntingTE(obj.Hz, obj.Ex, obj.Ey); 
                    
                case 'TM'
                    [Sx, Sy] = poyntingTM(obj.Ez, obj.Hx, obj.Hy);  
                    
            end
        end
        
        %% Set up the permittivity distribution 
        function obj = eps_setup(obj)     
            obj.dL = [diff(obj.xrange) diff(obj.yrange)]./obj.N; 
            
            obj.omega0 = 2*pi*obj.c0 / (obj.wvlen0 * obj.L0); 
            
            
            obj.eps_r = ones(obj.N); 

            N_blocks = length(obj.blocks); 
            for i = 1 : N_blocks
                obj.eps_r = assign_blocks(obj.eps_r, obj.blocks{i}{1}, obj.blocks{i}{2}, obj.blocks{i}{3}, obj.xrange, obj.yrange, obj.N); 
            end
            
        end
        
        %% Set up the source
        function [beta, obj] = source_setup(obj, src_type, src_loc, src_params)
        %% Output parameter: 
        %   beta: for a modal source, this outputs the propagation constant
        
        %%
            obj.src = zeros(obj.N); 
            beta = []; 
            
            src_type = lower(src_type); 
            switch src_type
                %% Point source
                case 'point'
                    % src_loc is the location of the point source
                    % src_params contains the amplitude of the source
                    src_ind_x = round((src_loc(1) - obj.xrange(1)) / diff(obj.xrange) * obj.N(1)) + 1; 
                    src_ind_y = round((src_loc(2) - obj.yrange(1)) / diff(obj.yrange) * obj.N(2)) + 1; 
                    
                    obj.src(src_ind_x, src_ind_y) = src_params; 
                    
                %% Total field scattered field
                case 'tfsf'
                    % src_loc contains [x_start, y_start, x_length, y_length]
                    % src_params contains [amplitude, angle] of the source
                    
                    hx = obj.dL(1); 
                    hy = obj.dL(2); 
                    
                    % Q = mask function. It is 1 in SF region, 0 in TF region
                    M = prod(obj.N); 
                    Q = ones(obj.N); 
                    within_TF = @(x, y) x > src_loc(1) & y > src_loc(2) & x < src_loc(1)+src_loc(3) & y < src_loc(2)+src_loc(4); 
                    Q = assign_val(Q, obj.xrange, obj.yrange, within_TF, 0); 
                    Q_sparse = spdiags(Q(:), 0, M, M); 
                    
                    % Calculate A matrix without the scatterer
                    switch obj.pol
                        case 'TM'
                            [~, ~, ~, A] = solveTM(obj.L0, obj.wvlen0, obj.xrange, obj.yrange, obj.eps_r, zeros(obj.N), obj.Npml); 
                        case 'TE'
                            [~, ~, ~, A] = solveTE(obj.L0, obj.wvlen0, obj.xrange, obj.yrange, obj.eps_r, zeros(obj.N), obj.Npml); 
                    end
                    
                    % Calculate the TFSF plane wave
                    theta = src_params(2) * pi/180; % Angle of incidence from X-AXIS
                    x_vec = linspace(obj.xrange(1)+hx/2, obj.xrange(2)-hx/2, obj.N(1)); 
                    y_vec = linspace(obj.yrange(1)+hy/2, obj.yrange(2)-hy/2, obj.N(2)); 

                    [X, Y] = meshgrid(x_vec, y_vec); 
                    f_src = src_params(1) * exp(-1i * 2*pi/obj.wvlen0 * (X'*cos(theta) + Y'*sin(theta))); 
                    f_src_vec = f_src(:)/obj.omega0; 
                    
                    % Calculate the TFSF source
                    b = (Q_sparse*A - A*Q_sparse) * f_src_vec; 
                    obj.src = reshape(b, obj.N); 
                    
                %% Waveguide modal source
                case 'modal'
                    % src_loc is the center coordinate of the modal source
                    % src_params is a cell that contains: 
                    %   {orientation, mode, num_points}
                    
                    src_ind_x = round((src_loc(1) - obj.xrange(1)) / diff(obj.xrange) * obj.N(1)) + 1; 
                    src_ind_y = round((src_loc(2) - obj.yrange(1)) / diff(obj.yrange) * obj.N(2)) + 1; 
                    
                    orientation = src_params{1}; 
                    mode = src_params{2}+1; 
                    n_pts = src_params{3}; 
                    
                    %% orientation = 'v' for vertical, 'h' for horizontal
                    orientation = lower(orientation); 
                    switch orientation
                        case 'v'
                            % Get the permittivity cross section
                            src_line_inds = src_ind_y-n_pts : src_ind_y+n_pts; 
                            eps_src = obj.eps_r(src_ind_x, src_line_inds); 
                            
                            src_xrange = [src_line_inds(1)*obj.dL(2), src_line_inds(end)*obj.dL(2)]; 
                            % Solve for the modes
                            eps_src = transpose(eps_src); 
                            [temp_src, beta] = solve_wg_modes_1D(obj.pol, obj.L0, obj.wvlen0, src_xrange, [-1 1], eps_src, mode); 
                                                        
                            % Assign the source profile
                            obj.src(src_ind_x, src_line_inds) = transpose(temp_src); 
                            
                        case 'h'
                            % Get the permittivity cross section
                            src_line_inds = src_ind_x-n_pts : src_ind_x+n_pts; 
                            eps_src = obj.eps_r(src_line_inds, src_ind_y); 
                            
                            src_xrange = [src_line_inds(1)*obj.dL(1), src_line_inds(end)*obj.dL(1)]; 
                            
                            % Solve for the modes
                            [temp_src, beta] = solve_wg_modes_1D(obj.pol, obj.L0, obj.wvlen0, src_xrange, [-1 1], eps_src, mode); 
                            
                            % Assign the source profile
                            obj.src(src_line_inds, src_ind_y) = temp_src; 
                    end
            end
        end
    end
    
end

