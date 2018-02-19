function Dws = createDws(w, s, dL, N)
%% This function creates the sparse derivative operators in FDFD
%% Input parameters
%   w: one of 'x', 'y'
%   s: one of 'f' and 'b'
%   dL: [dx dy] for 2D
%   N: [Nx Ny] for 2D
%
%% Output parameters
%   Dws: derivative operators

%% Compute Nx, Ny
Nx = N(1); 
Ny = N(2); 

%% Sparse identity matrices
Ix = speye(Nx); 
Iy = speye(Ny); 

%% Create derivative operators
switch w
    case 'x'
        if s == 'f'
            dxf = -Ix + circshift(Ix, [0 1]); 
            Dws = 1/dL(1) * kron(Iy, dxf); 
        else
            dxb = Ix - circshift(Ix, [0 -1]); 
            Dws = 1/dL(1) * kron(Iy, dxb); 
        end
        
        
    case 'y'
        if s == 'f'
            dyf = -Iy + circshift(Iy, [0 1]); 
            Dws = 1/dL(2) * kron(dyf, Ix); 
        else
            dyb = Iy - circshift(Iy, [0 -1]); 
            Dws = 1/dL(2) * kron(dyb, Ix); 
        end
end

end

