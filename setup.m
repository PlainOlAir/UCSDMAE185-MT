%% Air at Standard Conditions (sea-level)
p0 = 101300;    % Pressure (Pa)
T0 = 288.15;    % Temperature (K)
rho0 = 1.225;   % Density (kg/m^3)
R = 287;        % Gas constant (J/(kg-K))
cp = 1005;      % Isobaric specific heat capacity (J/(kg-K))
cv = 718;       % Volumetric specific heat capacity (J/(kg-K))
gamma = 1.4;    % Ratio of specific heats
mu0 = 1.735e-5;  % Dynamic viscosity (N-s/m^2)
S1 = 110.4;     % Sutherland's temperature (K)
Pr = 0.71;      % Prandtl number 

a = sqrt(gamma * R * T0);
uinf = a * M;

%% Grid setup
M = 4;      % Mach Number
nx = 75;    % x Number of points
ny = 80;    % y Number of points
L = 1e-5;   % Length of computational domain (m)
H = 8e-6;   % Height of computational domain (m)

dx = L / nx;
dy = H / ny;

step_total = 1500;

[x, y] = ndgrid(linspace(0, L, nx), linspace(0, H, ny));

%% Preallocation
[u,v,p,rho,T,e,U_bar] = deal(zeros(nx, ny));
U = zeros(4,nx,ny,step_total);

%% Initial Conditions
u(:, :) = uinf;
v(:, :) = 0;
p(:, :) = p0;
T(:, :) = T0;
rho(:,:) = rho0;

%% Boundary Conditions
% Bottom wall
u(:, 1) = 0;
v(:, 1) = 0;
T(:, 1) = 0;

% Inlet
u(1, :) = uinf;
v(1, :) = 0;
p(1, :) = p0;
T(1, :) = T0;

% Far field
u(:, end) = uinf;
v(:, end) = 0;
p(:, end) = p0;
T(:, end) = T0;

% in general:
% U = 4 x Nx x Ny x step array of conservative variables:
%   1 - density
%   2 - x-mass flux
%   3 - y-mass flux
%   4 - total energy
% set U(:,:,:,step 1) to BC's
U(:,:,:,1) = prim2cons(rho,u,v,T,cv);