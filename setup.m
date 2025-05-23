%% Air at Standard Conditions (sea-level)
p0 = 101300;                % Pressure (Pa)
T0 = 288.15;                % Temperature (K)
rho0 = 1.225;               % Density (kg/m^3)
R = 287;                    % Gas constant (J/(kg-K))
cp = 1005;                  % Isobaric specific heat capacity (J/(kg-K))
cv = 718;                   % Volumetric specific heat capacity (J/(kg-K))
gamma = 1.4;                % Ratio of specific heats
mu0 = 1.735e-5;             % Dynamic viscosity (N-s/m^2)
S1 = 110.4;                 % Sutherland's temperature (K)
Pr = 0.71;                  % Prandtl number 
M = 4;                      % Mach Number

a0 = sqrt(gamma * R * T0);  % Speed of Sound
uinf = a0 * M;              % Free-stream velocity
pinf = p0;
Tinf = T0;

%% Grid setup
nx = 75;    % x Number of points
ny = 80;    % y Number of points
L = 1e-5;   % Length of computational domain (m)
H = 8e-6;   % Height of computational domain (m)
dx = L / nx;
dy = H / ny;
step_total = 1500;
time = zeros(1,step_total);
[xx, yy] = ndgrid(linspace(0, L, nx), linspace(0, H, ny));

%% Preallocation
[u,v,p,rho,T,e,UBar] = deal(zeros(nx, ny));
% U = 4 x Nx x Ny x step array of conservative variables:
%   1 - density
%   2 - x-mass flux
%   3 - y-mass flux
%   4 - total energy
U = zeros(4,nx,ny);

%% Initial Conditions
u(:, :) = uinf;
v(:, :) = 0;
p(:, :) = pinf;
T(:, :) = Tinf;
[rho, u, v, T, p, e, Et] = bc_enforcer(u, v, T, p, cv, R, uinf, pinf, Tinf);
U(:,:,:) = prim2cons(rho,u,v,T,cv);
new_vars = {rho, u, v, e, p, T};
U_prev = U;
convergence = zeros(step_total+1, 4);
% output variables are rho, u, v, e, p, T, convergence
output_vars = cell(1,7);
for i = 1:6
    output_vars{i} = zeros(nx, ny, step_total+1);
    output_vars{i}(:, :, 1) = new_vars{i};
end
output_vars{7} = zeros(4,step_total);