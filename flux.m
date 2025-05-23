function [E, F] = flux(FD_method, U, dx, dy, mu, k, R, cv)
%% Description %%
% Calculates conservative variable fluxes in x and y directions as E and F
% respectively in the 2-D compressible Navier-Stokes equations
% 
% - DIRECTION OF FD MUST BE OPPOSITE OF THAT USED IN 'maccormack.m'

% INPUTS
% FD_method = finite difference method (backward or forward)
% U = conservative variable matrix with size of 4 x Nx x Ny
% dx = x-direction spacing of the data points
% dy = y-direction spacing of the data points
% mu = dynamic viscosity matrix with size Nx x Ny
% k = heat transfer coefficient matrix with size Nx x Ny
% R = specific gas constant of fluid
% cv = volumetric specific heat capacity of fluid

% OUTPUTS
% E = x-direction flux vector with size Nx x Ny
% F = y-direction flux vector with size Nx x Ny


%% Setup %%

% Preallocate E & F array
[~, nx, ny, ~] =  size(U);
[E, F] = deal(zeros(4, nx,ny));

%% Calculate %%

% Extract primitive variables
[rho, u, v, T, p, ~, Et] = cons2prim(U, R, cv);

% Apply forward or backwards finite difference
if strcmpi(FD_method,'forward')
    
    % Solve divergence of flow field
    div_u = ddx_fwd(u, dx, 1) + ddx_fwd(v, dy, 2);

    % Solve for shear stresses, E & F must have separate tau_xy to prevent
    % double biasing
    tau_xyE = mu .* (ddx_central(u, dy, 2) + ddx_fwd(v, dx, 1));
    tau_xyF = mu .* (ddx_fwd(u, dy, 2) + ddx_central(v, dx, 1));
    tau_xx = 2 * mu .* (ddx_fwd(u, dx, 1) - (1/3) .* div_u);
    tau_yy = 2 * mu .* (ddx_fwd(v, dy, 2) - (1/3) .* div_u);

    % Solve for heat transfer rates
    qx = -k .* ddx_fwd(T, dx, 1);
    qy = -k .* ddx_fwd(T, dy, 2);

elseif strcmpi(FD_method,'backward')

    % Solve divergence of flow field
    div_u = ddx_bwd(u, dx, 1) + ddx_bwd(v, dy, 2);

    % Solve for shear stresses, E & F must have separate tau_xy to prevent
    % double biasing
    tau_xyE = mu .* (ddx_central(u, dy, 2) + ddx_bwd(v, dx, 1));
    tau_xyF = mu .* (ddx_bwd(u, dy, 2) + ddx_central(v, dx, 1));
    tau_xx = 2 .* mu .* (ddx_bwd(u, dx, 1) - (1/3) .* div_u);
    tau_yy = 2 .* mu .* (ddx_bwd(v, dy, 2) - (1/3) .* div_u);

    % Solve for heat transfer rates
    qx = -k .* ddx_bwd(T, dx, 1);
    qy = -k .* ddx_bwd(T, dy, 2);
else
    error("Invalid FD method given. Must be 'forward' or 'backward'. ")
end

% Build E & F arrays
E(1,:,:) = rho .* u;
E(2,:,:) = rho .* u .^ 2 + p - tau_xx;
E(3,:,:) = rho .* u .* v - tau_xyE;
E(4,:,:) = (Et + p) .* u - u .* tau_xx - v .* tau_xyE + qx;

F(1,:,:) = rho .* v;
F(2,:,:) = rho .* u .* v - tau_xyF;
F(3,:,:) = rho .* v .^ 2 + p - tau_yy;
F(4,:,:) = (Et + p) .* v - v .* tau_yy - u .* tau_xyF + qy;

end