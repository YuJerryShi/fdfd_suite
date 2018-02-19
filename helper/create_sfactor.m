function sfactor_array = create_sfactor(wrange, s, omega, eps0, mu0, Nw, Nw_pml)
%% Input Parameters
% wrange: [wmin wmax], range of domain in w-direction including PML
% s: 'b' or 'f', indicating whether s-factor is for Dwb or Dwf
% omega: angular frequency
% eps0: vacuum permittivity
% mu0: vacuum permeability
% Nw: number of cells in w-direction
% Nw_pml: number of cells in PML

%% Output Parameter
% sfactor_array: 1D array with Nw elements containing PML s-factors for Dws

eta0 = sqrt(mu0/eps0);  % vacuum impedance
m = 3.5;  % degree of polynomial grading
lnR = -12;  % R: target reflection coefficient for normal incidence

% find dw 
hw = diff(wrange) / (Nw); 
dw = Nw_pml * hw; 

% Sigma function
sig_max = -(m+1) * lnR / (2*eta0*dw); 
sig_w = @(l) sig_max *(l/dw)^m; 

S = @(l) 1 - 1i*sig_w(l) / (omega*eps0); 

%% PML vector
% initialize sfactor_array
sfactor_array = ones(Nw, 1); 
for i = 1:Nw
    switch s
        case 'f'
            if i <= Nw_pml
                sfactor_array(i) = S(hw * (Nw_pml - i + 0.5)); 
                
            elseif i > Nw - Nw_pml
                sfactor_array(i) = S(hw * (i - (Nw - Nw_pml) - 0.5)); 
            end
                
        case 'b'
            if i <= Nw_pml
                sfactor_array(i) = S(hw * (Nw_pml - i + 1)); 
                
            elseif i > Nw - Nw_pml
                sfactor_array(i) = S(hw * (i - (Nw - Nw_pml) - 1)); 
            end
    end

end

