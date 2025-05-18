%% Air at Standard Conditions (sea-level)
p0 = 101300;    % Pressure (Pa)
T0 = 288.15;    % Temperature (K)
rho0 = 1.225;   % Density (kg/m^3)
R = 287;        % Gas constant (J/(kg-K))
cp = 1005;      % Isobaric specific heat capacity (J/(kg-K))
cv = 718;       % Volumetric specific heat capacity (J/(kg-K))
gamma = 1.4;    % Ratio of specific heats
u0 = 1.735e-5;  % Dynamic viscosity (N-s/m^2)
S1 = 110.4;     % Sutherland's temperature (K)
Pr = 0.71;      % Prandtl number 

%% Grid setup
M = 4;      % Mach Number
nx = 75;    % x Number of points
ny = 80;    % y Number of points
L = 1e-5;   % Length of computational domain (m)
H = 8e-6;   % Height of computational domain (m)

dx = L / nx;
dy = H / ny;

[x, y] = ndgrid(linspace(0, L, nx), linspace(0, H, ny));
u = zeros(nx, ny);
v = u;
p = u;
rho = u;
T = u;
e = u;

a = sqrt(gamma * R * T0);
uinf = a * M;

step_total = 1500;

%% Initial Conditions
u(:, :) = uinf;
v(:, :) = 0;
p(:, :) = p0;
T(:, :) = T0;

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

%% Preallocation
[U, Ubar, Upred] = deal(zeros(nx,ny,step_total));
