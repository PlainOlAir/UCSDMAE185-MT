function [rho, u, v, T, p, e, Et] = bc_enforcer(u, v, T, p, cv, R, uinf, pinf, Tinf)
%% Description %%
% Enforces boundary conditions on the computational domain
%
% - SHARP LEADING EDGE
% - ENFORCES INLET BOUNDARIES ON THE TOP FAR-FIELD (ENSURE SHOCK DOES NOT
%   TOUCH UPPER DOMAIN
% - 2nd ORDER ACCURATE EXTRAPOLATION AT OUTLET

% INPUTS
% u = x-velocity matrix with size Nx x Ny
% v = y-velocity matrix with size Nx x Ny
% T = temperature matrix with size Nx x Ny
% p = pressure matrix with size Nx x Ny
% cv = volumetric specific heat capacity of fluid
% R = specific gas constant
% uinf = inlet fluid velocity
% pinf = inlet fluid pressure
% Tinf = inlet fluid temperature

% OUTPUTS
% rho = density recalculated with BCs applied with size Nx x Ny
% u = x-velocity recalculated with BCs applied with size Nx x Ny
% v = y-velocity recalculated with BCs applied with size Nx x Ny
% T = temperature recalculated with BCs applied with size Nx x Ny
% p = pressure recalculated with BCs applied with size Nx x Ny
% e = internal energy recalculated with BCs applied with size Nx x Ny
% Et = total energy recalculated with BCs applied with size Nx x Ny
 
%% Process BCs %%

% Bottom wall
u(2:end, 1) = 0;
v(2:end, 1) = 0;
T(2:end, 1) = Tinf;
p(2:end, 1) = 2*(p(2:end,2)) - p(2:end,3);

% Inlet
u(1, 2:end) = uinf;
v(1, 2:end) = 0;
p(1, 2:end) = pinf;
T(1, 2:end) = Tinf;

% Far field
u(:, end) = uinf;
v(:, end) = 0;
p(:, end) = pinf;
T(:, end) = Tinf;

% Outlet
u(end,:) = 2*(u(end-1,:)) - u(end-2,:);
v(end,:) = 2*(v(end-1,:)) - v(end-2,:);
p(end,:) = 2*(p(end-1,:)) - p(end-2,:);
T(end,:) = 2*(T(end-1,:)) - T(end-2,:);

% Leading edge
u(1,1) = 0;
v(1,1) = 0;
p(1,1) = pinf;
T(1,1) = Tinf;

rho = p./(R.*T);                    % Density
e = cv.*T;                          % Internal energy
Et = rho .* (e + (u.^2 + v.^2)./2); % Total energy

end