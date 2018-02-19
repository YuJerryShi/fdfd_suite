classdef fdfd_blochX_modes < fdfd_main
    %FDFD_BLOCHX_SOLVE: Specializes in 1D periodic photonic crystal modes
    
    properties
        Kx % Bloch wave vector in x
        n_modes % Number of modes solved
        omega_eigs % Eigen-frequencies of the Bloch modes
        lambda_eigs % The corresponding wavelengths of the Eigen-frequencies
    end
    
    methods
        %% Run the simulation
        function obj = simulate(obj)
            switch obj.pol
                case 'TE'
                    [obj.Hz, obj.Ex, obj.Ey, obj.omega_eigs] = solveTE_BlochX(obj.L0, obj.wvlen0, obj.xrange, obj.yrange, obj.eps_r, obj.Kx, obj.Npml, obj.n_modes); 
                case 'TM'
                    [obj.Ez, obj.Hx, obj.Hy, obj.omega_eigs] = solveTM_BlochX(obj.L0, obj.wvlen0, obj.xrange, obj.yrange, obj.eps_r, obj.Kx, obj.Npml, obj.n_modes); 
            end
            
            obj.lambda_eigs = 2*pi*obj.c0/obj.L0 ./ obj.omega_eigs; 
        end
        
        %% Filter out unwanted modes of the simulation domain
        function obj = filtering(obj, y_planes, thresh)
            % thresh is a value between 0 and 1 
            
            % Get the indices of the box
            y_ind1 = floor((y_planes(1)-obj.yrange(1))/obj.dL(2)) + 1; 
            y_ind2 = floor((y_planes(2)-obj.yrange(1))/obj.dL(2)) + 1;  

            % Filtering to keep the desired modes
            omega_eigs_temp = []; 
            
            switch obj.pol
                case 'TE'
                    Hz_temp = {}; 
                    Ex_temp = {}; 
                    Ey_temp = {}; 
                    
                    % Apply the filter for each Hz
                    counter = 1; 
                    for i = 1:length(obj.Hz)
                        intensity1 = sum(sum(abs(obj.Hz{i}(:, y_ind1:y_ind2)))); 
                        intensity2 = sum(sum(abs(obj.Hz{i}))); 
                        
                        if intensity1 > thresh*intensity2
                            Hz_temp{counter, 1} = obj.Hz{i}; 
                            Ex_temp{counter, 1} = obj.Ex{i}; 
                            Ey_temp{counter, 1} = obj.Ey{i}; 
                            
                            omega_eigs_temp = [omega_eigs_temp; obj.omega_eigs(i)]; 
                            counter = counter + 1; 
                        end
                        
                    end
                    
                    % Override Hz, Ex, Ey, omega_eigs after filtering
                    obj.Hz = Hz_temp; 
                    obj.Ex = Ex_temp; 
                    obj.Ey = Ey_temp; 
                    obj.omega_eigs = omega_eigs_temp; 
                    
                case 'TM'
                    Ez_temp = {}; 
                    Hx_temp = {}; 
                    Hy_temp = {}; 
                                        
                    % Apply the filter for each Ez
                    counter = 1; 
                    for i = 1:length(obj.Ez)
                        intensity1 = sum(sum(abs(obj.Ez{i}(:, y_ind1:y_ind2)))); 
                        intensity2 = sum(sum(abs(obj.Ez{i}))); 
                        
                        if intensity1 > thresh*intensity2
                            Ez_temp{counter, 1} = obj.Ez{i};  
                            Hx_temp{counter, 1} = obj.Hx{i}; 
                            Hy_temp{counter, 1} = obj.Hy{i}; 
                            
                            omega_eigs_temp = [omega_eigs_temp; obj.omega_eigs(i)]; 
                            counter = counter + 1; 
                        end
                    end
                    
                    % Override Ez, Hx, Hy, omega_eigs after filtering
                    obj.Ez = Ez_temp; 
                    obj.Hx = Hx_temp; 
                    obj.Hy = Hy_temp; 
                    obj.omega_eigs = omega_eigs_temp; 
            end
            
            obj.lambda_eigs = 2*pi*obj.c0 ./ obj.omega_eigs; 
        end
        
    end
    
end

