function U = prim2cons(rho, u, v, T, cv)
%% Description %%
% Converts primitive variables into conservative variables
% Assumes calorically perfect gas

% INPUTS
% rho = density with size Nx x Ny
% u = x-velocity with size Nx x Ny
% v = y-velocity with size Nx x Ny
% T = temperature with size Nx x Ny
% cv = volumetric specific heat capacity

% OUTPUTS
% U = 4 x Nx x Ny array of conservative variables:
%   1 - density
%   2 - x-mass flux
%   3 - y-mass flux
%   4 - total energy

%% Setup %%
    [nx, ny] = size(rho); % get size for U array
    U = zeros(4, nx, ny); % preallocate memory

%% Calculate %%
    e = cv .* T;                        % internal energy
    Et = rho .* (e + (u.^2 + v.^2)./2); % total energy
    
    % build U array
    U(1, :, :) = rho;
    U(2, :, :) = rho .* u;
    U(3, :, :) = rho .* v;
    U(4, :, :) = Et;

end