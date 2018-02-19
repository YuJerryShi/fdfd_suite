classdef fdfd_wg_modes < fdfd_main
    %FDFD_WG_MODES: Specializes in solving for the guided modes of a
    %optical waveguide given its x-y cross section
    
    properties
        beta % Propagation constants for each waveguide mode in 1/L0
        n_modes % Number of waveguide modes
    end
    
    methods
        %% Run the simulation
        function obj = simulate(obj)
            [obj.Hx, obj.Hy, obj.Hz, obj.Ex, obj.Ey, obj.Ez, obj.beta] = solve_wg_modes(obj.L0, obj.wvlen0, obj.xrange, obj.yrange, obj.eps_r, obj.n_modes); 
        end
        
        %% Visualize the real part of the fields
        function [] = visreal_wg(obj, mode)
            subplot(2, 3, 1); 
            vis_real(obj.Hx{mode}, obj.xrange, obj.yrange); 
            title('Hx'); 

            subplot(2, 3, 2); 
            vis_real(obj.Hy{mode}, obj.xrange, obj.yrange); 
            title('Hy'); 

            subplot(2, 3, 3); 
            vis_real(obj.Hz{mode}, obj.xrange, obj.yrange); 
            title('Hz'); 

            subplot(2, 3, 4); 
            vis_real(obj.Ex{mode}, obj.xrange, obj.yrange); 
            title('Ex'); 

            subplot(2, 3, 5); 
            vis_real(obj.Ey{mode}, obj.xrange, obj.yrange); 
            title('Ey'); 

            subplot(2, 3, 6); 
            vis_real(obj.Ez{mode}, obj.xrange, obj.yrange); 
            title('Ez'); 
        end
        
        %% Visualize the amplitude of the fields
        function [] = visabs_wg(obj, mode)
            subplot(2, 3, 1); 
            vis_abs(obj.Hx{mode}, obj.xrange, obj.yrange); 
            title('Hx'); 

            subplot(2, 3, 2); 
            vis_abs(obj.Hy{mode}, obj.xrange, obj.yrange); 
            title('Hy'); 

            subplot(2, 3, 3); 
            vis_abs(obj.Hz{mode}, obj.xrange, obj.yrange); 
            title('Hz'); 

            subplot(2, 3, 4); 
            vis_abs(obj.Ex{mode}, obj.xrange, obj.yrange); 
            title('Ex'); 

            subplot(2, 3, 5); 
            vis_abs(obj.Ey{mode}, obj.xrange, obj.yrange); 
            title('Ey'); 

            subplot(2, 3, 6); 
            vis_abs(obj.Ez{mode}, obj.xrange, obj.yrange); 
            title('Ez'); 
        end
        
        %% Draw quiver plots
        function [] = quiver_wg(obj, mode)
            x_total = 1:obj.N(1); 
            y_total = 1:obj.N(2); 

            resolution = 6; 
            indx = round(resolution/2):resolution:obj.N(1); 
            indy = round(resolution/2):resolution:obj.N(2); 

            x = x_total(indx); 
            y = y_total(indy); 

            [X_temp, Y_temp] = meshgrid(x, y); 
            X = X_temp/obj.N(1) * diff(obj.xrange) + obj.xrange(1); 
            Y = Y_temp/obj.N(2) * diff(obj.yrange) + obj.yrange(1); 

            Ex_qui = real(obj.Ex{mode}(indx, indy)); 
            Ey_qui = real(obj.Ey{mode}(indx, indy)); 

            Hx_qui = real(obj.Hx{mode}(indx, indy)); 
            Hy_qui = real(obj.Hy{mode}(indx, indy)); 
            
            subplot(1, 2, 1); 
            quiver(X, Y, transpose(Ex_qui), transpose(Ey_qui), 'AlignVertexCenters', 'on');
            axis equal
            axis([obj.xrange obj.yrange]); 
            title('E field'); 

            subplot(1, 2, 2); 
            
            quiver(X, Y, transpose(Hx_qui), transpose(Hy_qui), 'AlignVertexCenters', 'on'); 
            
            axis equal
            axis([obj.xrange obj.yrange]); 
            title('H field'); 
        end
        
    end
    
end

