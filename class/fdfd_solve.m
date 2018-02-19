classdef fdfd_solve < fdfd_main
    %FDFD_SOLVE: specializes in a direct solution for a photonic structure
    %under the excitation of a source
    
    properties
        A % Sparse system matrix of the solver
    end
    
    methods
        %% Run the simulation
        function obj = simulate(obj)
            switch obj.pol
                case 'TE'
                    [obj.Hz, obj.Ex, obj.Ey, obj.A] = solveTE(obj.L0, obj.wvlen0, obj.xrange, obj.yrange, obj.eps_r, obj.src, obj.Npml); 
                case 'TM'
                    [obj.Ez, obj.Hx, obj.Hy, obj.A] = solveTM(obj.L0, obj.wvlen0, obj.xrange, obj.yrange, obj.eps_r, obj.src, obj.Npml); 
                    
            end
        end
        
    end
    
end

