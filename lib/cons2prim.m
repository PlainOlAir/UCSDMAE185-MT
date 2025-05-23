function [rho, u, v, T, p, e, Et] = cons2prim(U, R, cv, nx, ny)
%% Description %%
% Converts conservative variables into primitive variables
% Assumes calorically perfect gas

% INPUTS
% U = conservative varibles with size 4 x Nx x Ny
% R = specific gas constant
% cv = volumetric specific heat capacity
% nx = number of rows in U (OPTIONAL BUT RECCOMENDED)
% ny = number of columns in U (OPTIONAL BUT RECCOMENDED)

% OUTPUTS
% rho = density with size Nx x Ny
% u = x-velocity with size Nx x Ny
% v = y-velocity with size Nx x Ny
% T = temperature with size Nx x Ny
% p = pressure with size Nx x Ny
% e = internal energy with size Nx x Ny
% Et = total energy with size Nx x Ny

%% Setup %%

% Process input variables
if nargin == 3
    [~, nx, ny] = size(U);
end

%% Calculate %%

% Separate values from conservative variable array

rho = reshape(U(1, :, :), [nx, ny]);        % Density
u   = reshape(U(2, :, :), [nx, ny]) ./ rho; % X-velocity
v   = reshape(U(3, :, :), [nx, ny]) ./ rho; % Y-velocity
Et  = reshape(U(4, :, :), [nx, ny]);        % Total energy

e   = Et ./ rho - (u.^2 + v.^2) ./ 2;       % Internal energy

T   = e ./ cv;                              % Temperature

p   = rho .* R .* T;                        % Pressure

end