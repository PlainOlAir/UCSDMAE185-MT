clear

addpath('lib\')

setup

dt = 2.35e-11; % Constant time step (s)

time = 0;

for step = 1:step_total
    % primitive var's and U update
    U(:,:,step) = prim2cons(rho(:, :, step), u(:, :, step), v(:, :, step), T(:, :, step), cv);
    % I/O and visualization
    % compute delta_t CFL
    %% Predictor
    if mod(step, 2) == 1
        FD_method = 'fwd';
    else
        FD_method = 'bwd';
    end
    [E(:,:,step), F(:,:,step)] = flux(U(:,:,step),dx,dy,FD_method);
    U_bar = U(:,:,step) - dt.*(E(:,:,step) + F(:,:,step));
    [rho_bar, u_bar, v_bar, T_bar, p_bar, e_bar, Et_bar] = cons2prim(U_bar, R, cv);
    
    % Update U using the predicted values
    U_pred = prim2cons(rho_bar, u_bar, v_bar, T_bar, cv);
    
    %% Corrector
    if mod(step, 2) == 1

    else
    end
end

function [E, F] = flux(U, dx, dy, FD_method) 
% need gradient and bias eq's
    [rho, u, v, T, p, ~, Et] = cons2prim(U, R, cv);
    if strcmpi(FD_method,'fwd') == 1 % predictor
        tau_xy = mu(ddx_fwd(u,dy,2,true) + ddx_fwd(v,dx,1,true));
        tau_xx = 2.*(mu.*ddx_fwd(u,dx,1,true) - (1/3).*gradient.*u);
        tau_yy = 2.*(mu.*ddx_fwd(v,dy,2,true) - (1/3).*gradient.*v);
        qx = -k.*ddx_fwd(T,dx,1,true);
        qy = -k.*ddx_fwd(T,dy,2,true);
    elseif strcmpi(FD_method,'bwd') == 1 % corrector
        tau_xy = mu(ddx_bwd(u,dy,2,true) + ddx_bwd(v,dx,1,true));
        tau_xx = 2.*(mu.*ddx_bwd(u,dx,1,true) - (1/3).*gradient.*u);
        tau_yy = 2.*(mu.*ddx_bwd(v,dy,2,true) - (1/3).*gradient.*v);
        qx = -k.*ddx_bwd(T,dx,1,true);
        qy = -k.*ddx_bwd(T,dy,2,true);
    elseif strcmpi(FD_method,'fwd_bias') == 1
        % tau_xy = mu(ddx_bwd(u,dy,2,true) + ddx_bwd(v,dx,1,true));
        % tau_xx = 2.*(mu.*ddx_bwd(u,dx,1,true) - (1/3).*gradient.*u);
        % tau_yy = 2.*(mu.*ddx_bwd(v,dy,2,true) - (1/3).*gradient.*v);
        % qx = -k.*ddx_bwd(T,dx,1,true);
        % qy = -k.*ddx_bwd(T,dy,2,true);
    elseif strcmpi(FD_method,'bwd_bias') == 1
        % tau_xy = mu(ddx_bwd(u,dy,2,true) + ddx_bwd(v,dx,1,true));
        % tau_xx = 2.*(mu.*ddx_bwd(u,dx,1,true) - (1/3).*gradient.*u);
        % tau_yy = 2.*(mu.*ddx_bwd(v,dy,2,true) - (1/3).*gradient.*v);
        % qx = -k.*ddx_bwd(T,dx,1,true);
        % qy = -k.*ddx_bwd(T,dy,2,true);
    else
        error('Incorrect or no FD method given.')
    end
    E(1,1) = rho.*u;
    E(2,1) = rho.*u.^2 + p - tau_xx;
    E(3,1) = rho.*u.*v - tau_xy;
    E(4,1) = (Et + p).*u - u.*tau_xx - v.*tau_xy + qx;
    F(1,1) = rho.*v;
    F(2,1) = rho.*u.*v - tau_xy;
    F(3,1) = rho.*v.^2 + p - tau_yy;
    F(4,1) = (Et + p).*v - v.*tau_yy - u.*tau_xy + qy;
end