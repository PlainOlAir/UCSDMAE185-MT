function dx3 = dxn_3d(func, f, dx, dir, periodicBC)
%% Description %%
% Finite difference functions cannot handle 3-D arrays in the form of 
% Nvars x Nx x Ny
% which is common when dealing with conservative variables
% This function loops through Nvars and applies the finite differences to
% each 2-D array of the variables

% INPUTS
% func = function to apply to f (ddx_bwd, ddx_fwd, ddx_central, d2dx2)
% f = data matrix with the size of Nvars x Nx x Ny
% dx = x-direction spacing of the data points
% dir = direction (1 = x, 2 = y) (OPTIONAL)
% periodicBC = true/false to set periodic boundaries (OPTIONAL)

% OUTPUTS
% dx3 = derivative "func" applied to each variable in f

%% Setup %%

% Process input variables
switch nargin
    case 3
        dir = 1;
        periodicBC = false;
    case 4
        periodicBC = false;
end

% Preallocate output array
[nvars, nx, ny] = size(f);
dx3 = zeros(nvars, nx, ny);

%% Calculate

% Loop through variables and apply func
for i = 1:nvars
    dx3(i, :, :) = func(reshape(f(i, :, :), [nx, ny]), dx, dir, periodicBC);
end

end