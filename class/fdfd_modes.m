classdef fdfd_modes < fdfd_main
    %FDFD_MODES: Specializes in solving for the eigenmodes of 2D
    %photonic structures
    
    properties
        n_modes % Number of modes 
        omega_eigs % Eigen-frequencies of the modes
        lambda_eigs % Wavelengths corresponding to the eigen-frequencies
        A % Sparse system matrix used for eigenmode solution
    end
    
    methods
        %% Run the simulation
        function obj = simulate(obj)
            switch obj.pol
                case 'TE'
                    [obj.Hz, obj.Ex, obj.Ey, obj.A, obj.omega_eigs] = solveTE_modes(obj.L0, obj.wvlen0, obj.xrange, obj.yrange, obj.eps_r, obj.Npml, obj.n_modes); 
                case 'TM'
                    [obj.Ez, obj.Hx, obj.Hy, obj.A, obj.omega_eigs] = solveTM_modes(obj.L0, obj.wvlen0, obj.xrange, obj.yrange, obj.eps_r, obj.Npml, obj.n_modes); 
            end
            
            obj.lambda_eigs = 2*pi*obj.c0/obj.L0 ./ obj.omega_eigs; 
        end
        
    end
    
end

