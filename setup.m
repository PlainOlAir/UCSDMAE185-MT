%%%%%%%%%%%%%%%%%%%%%%
%%% Setup Function %%%
%%%%%%%%%%%%%%%%%%%%%%

%% Data & Animation
step_total = 1500;                  % Maximum number of simulation steps

animationstep = 10;                 % Steps between frames
animationwidth = 1280;              % Animation window width (pixels)
animationheight = 720;              % Animation window height (pixels)
animationexport = false;             % Set to true if export to gif
animationfile = 'animation.gif';    % Animiation gif file name
animationtime = 5;                  % Animation time (s)

%% Grid setup
nx = 75;    % x Number of points
ny = 80;    % y Number of points
L = 1e-5;   % Length of computational domain (m)
H = 8e-6;   % Height of computational domain (m)

dx = L / nx; % Find x cell spacing (m)
dy = H / ny; % Find y cell spacing (m)

% Build x & y coordinate grid arrays
[xx, yy] = ndgrid(linspace(0, L, nx), linspace(0, H, ny));

%% Preallocation
% Build primitive arrays to enforce initial and boundary conditions at
% first time step
[rho, u, v, T, p, e, Et] = deal(zeros(nx, ny)); 

% Preallocate output data array & time array
output_vars = cell(1,7);
[output_vars{1:6}] = deal(zeros(nx, ny, step_total+1));
output_vars{7} = zeros(step_total+1, 4);
time = zeros(1,step_total+1);

%% Air Conditions (Default: Standard Conditions (sea-level))
p0      = 101300;   % Pressure (Pa)
T0      = 288.15;   % Temperature (K)
rho0    = 1.225;    % Density (kg/m^3)
R       = 287;      % Gas constant (J/(kg-K))
cp      = 1005;     % Isobaric specific heat capacity (J/(kg-K))
cv      = 718;      % Volumetric specific heat capacity (J/(kg-K))
gamma   = 1.4;      % Ratio of specific heats
mu0     = 1.735e-5; % Dynamic viscosity (N-s/m^2)
S1      = 110.4;    % Sutherland's temperature (K)
Pr      = 0.71;     % Prandtl number 
M       = 4;        % Mach Number

%% Initial Conditions
time = 0;                   % Time (s)
a0 = sqrt(gamma * R * T0);  % Speed of Sound (m/s)
uinf = a0 * M;              % Free-stream velocity (m/s)
pinf = p0;                  % Free-stream pressure (Pa)
Tinf = T0;                  % Free-stream temperature (K)

% Initialize entire flow field
u(:, :) = uinf; 
v(:, :) = 0;
p(:, :) = pinf;
T(:, :) = Tinf;

%% Boundary Conditions
% Enforce BCs with 'bc_enforcer.m'
[rho, u, v, T, p, e, Et] = bc_enforcer(u, v, T, p, cv, R, uinf, pinf, Tinf);

%% Final Setup
% Build initial U with ICs and BCs
U = prim2cons(rho,u,v,T,cv);
new_vars = {rho, u, v, e, p, T};
U_prev = U;

% Store initial values
output_vars = cell(1,6);
for i = 1:6
    output_vars{i} = zeros(nx, ny, step_total+1);
    output_vars{i}(:, :, 1) = new_vars{i};
end