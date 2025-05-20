function [rho, u, v, T, p, e, Et] = cons2prim(U, R, cv)
%% Description
% Converts conservative variables into primitive variables
% Assumes calorically perfect gas

% INPUTS
% U = conservative varibles with size 4 x Nx x Ny
% R = specific gas constant
% cv = volumetric specific heat capacity

% OUTPUTS
% rho = density with size Nx x Ny
% u = x-velocity with size Nx x Ny
% v = y-velocity with size Nx x Ny
% T = temperature with size Nx x Ny
% p = pressure with size Nx x Ny
% e = internal energy with size Nx x Ny
% Et = total energy with size Nx x Ny

%% Calculate
    % Separate values from conservative variable array
    rho = squeeze(U(1, :, :));
    u = squeeze(U(2, :, :)) ./ rho;
    v = squeeze(U(3, :, :)) ./ rho;
    Et = squeeze(U(4, :, :));

    e = squeeze(Et ./ rho - (u.^2 + v.^2) ./ 2); % internal energy

    T = squeeze(e ./ cv); % temperature

    p = squeeze(rho .* R .* T); % pressure

end