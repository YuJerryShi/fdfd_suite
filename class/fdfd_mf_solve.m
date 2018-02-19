classdef fdfd_mf_solve < fdfd_main
    %FDFD_MODE_SOLVE: Specializes in solving field profiles in
    %time-modulated photonic structures
    % Details of the implementation: Y. Shi, et al. Optica 3, pp.1256 (2016)
    
    properties
        Nsb % Number of sideband frequencies
        mod_amp_blocks % Cells containing modulation amplitude geometries
        mod_phi_blocks % Cells containing modulation phase geometries
        mod_amp % Modulation amplitude distribution
        mod_phi % Modulation phase distribution
        Omega % Modulation frequency
        omega % List of sideband frequencies
    end
    
    methods
        %% Visualize the modulation amplitude and phase
        function obj = vis_mod(obj)
            subplot(2, 1, 1); 
            vis_abs(obj.mod_amp, obj.xrange, obj.yrange); 
            title('Mod amplitude'); 
            
            subplot(2, 1, 2); 
            vis_real(obj.mod_phi, obj.xrange, obj.yrange); 
            title('Mod phase'); 
        end
        
        %% Run the simulation
        function obj = simulate(obj)
            switch obj.pol
                case 'TE'
                    [obj.Hz, obj.Ex, obj.Ey, obj.omega] = solveTE_mf(obj.L0, obj.wvlen0, obj.Omega, obj.Nsb, obj.xrange, obj.yrange, obj.eps_r, obj.mod_amp, obj.mod_phi, obj.src, obj.Npml); 
                case 'TM'
                    [obj.Ez, obj.Hx, obj.Hy, obj.omega] = solveTM_mf(obj.L0, obj.wvlen0, obj.Omega, obj.Nsb, obj.xrange, obj.yrange, obj.eps_r, obj.mod_amp, obj.mod_phi, obj.src, obj.Npml); 
            end
        end
        
        
        %% Add a modulation amplitude region
        function obj = add_mod_amp(obj, block_name, block_loc, block_val)
            N_blocks = length(obj.mod_amp_blocks); 

            obj.mod_amp_blocks{N_blocks+1} = {block_name, block_loc, block_val}; 
        end
        
        %% Add a modulation phase region
        function obj = add_mod_phi(obj, block_name, block_loc, block_val)
            N_blocks = length(obj.mod_phi_blocks); 

            obj.mod_phi_blocks{N_blocks+1} = {block_name, block_loc, block_val}; 
        end
        
        %% Set up the modulation region 
        function obj = mod_setup(obj)     
            obj.mod_amp = zeros(obj.N); 
            obj.mod_phi = zeros(obj.N); 
            
            % Assign the modulation amplitude matrix
            N_amp = length(obj.mod_amp_blocks); 
            for i = 1 : N_amp
                obj.mod_amp = assign_blocks(obj.mod_amp, obj.mod_amp_blocks{i}{1}, obj.mod_amp_blocks{i}{2}, obj.mod_amp_blocks{i}{3}, obj.xrange, obj.yrange); 
            end
            
            % Assign the modulation phase matrix
            N_phi = length(obj.mod_phi_blocks); 
            for i = 1 : N_phi
                obj.mod_phi = assign_blocks(obj.mod_phi, obj.mod_phi_blocks{i}{1}, obj.mod_phi_blocks{i}{2}, obj.mod_phi_blocks{i}{3}, obj.xrange, obj.yrange); 
            end
        end
        
    end
    
end

