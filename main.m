clear
close all
clc

addpath('lib\')

setup;

dt = 2.35e-11; % Constant time step (s)

time = 0;

%% --- Main Loop --- TODO: add BCs, I/O
for step = 1:step_total
    %% I/O, loop updates, delta_t_CFL, visualization
    a = sqrt(gamma*R*T);
    % rho, u, v, e, p, T, convergence
    convergence = p_previous - p;
    output_vars(:, :, :, step) = permute(cat(3, rho, u, v, e, p, T, convergence), [3 1 2]);
    p_previous = p;
    % compute delta_t CFL
    time(step+1,1) = time(step,1) + dt;
    %% Alternator
    if mod(step, 2) == 1
        FD_method_pred = 'fwd';
        FD_method_corr = 'bwd';
    else
        FD_method_pred = 'bwd_bias';
        FD_method_corr = 'fwd_bias';
    end
    %% Predictor
    % update mu, k
    mu = sutherland(T, mu0, T0, S1);
    k = (cp/Pr)*mu;
    % compute derivatives, update E, F
    [E, F] = flux(U(:,:,:,step),dx,dy,FD_method_pred,R,cv,mu,k);
    % compute U_bar using U, E, F
    U_bar = U(:,:,:,step) - dt*(E + F);
    % compute primitive var's from U_bar
    [rho_bar, u_bar, v_bar, T_bar, p_bar, e_bar, Et_bar] = cons2prim(U_bar, R, cv); 
    % enforce BC's on p, u,v, T (update rho, e,...)
    [p, u, v, T, rho, e, Et] = bc_enforcer(p_bar, u_bar, v_bar, T_bar, cv, R, uinf, p0, T0);
    U_bar = prim2cons(rho, u, v, T, cv);

    %% Corrector
    % update mu, k
    mu = sutherland(T, mu0, T0, S1);
    k = (cp/Pr)*mu;
    % compute derivatives, update E, F
    [E, F] = flux(U_bar,dx,dy,FD_method_corr, R, cv, mu, k);
    % compute U from primitive vars
    U(:,:,:,step+1) = (0.5).*(U(:,:,:,step) + U_bar - dt.*(E + F));
    % compute primitive vars from U
    [rho, u, v, T, p, e, Et] = cons2prim(U(:,:,:,step), R, cv);
    % enforce BC's on p, u, v, T (update rho, e,...)
    [p, u, v, T, rho, e, Et] = bc_enforcer(p, u, v, T, cv, R, uinf, p0, T0);
    % compute U from primitive vars
    U(:,:,:,step+1) = prim2cons(rho,u,v,T,cv);
 end
%% --- Animator --- TODO: verify functioning
figure;
tile = tiledlayout(2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
% initialize handles
axesArray = gobjects(1, 8);
h = gobjects(1, 8);
var_labels = {'$\rho$', '$u$', '$v$', '$e$', '$p$', '$T$', '$Convergence$'};
titles = gobjects(1, 8);
% initial plot for each variable
for var = 1:7
    axesArray(var) = nexttile(tile, var);
    h(var) = pcolor(axesArray(var), xx, yy, squeeze(output_vars(var,:,:,1)));
    shading(axesArray(var), 'interp');
    axis(axesArray(var), 'equal', 'tight');
    xlabel(axesArray(var), '$x$', 'Interpreter','latex');
    ylabel(axesArray(var), '$y$', 'Interpreter','latex');
    % colorbar(axesArray(var));
    titles(var) = title(axesArray(var), ...
        sprintf('%s at $t=%.3f$ s \\(%d/%d\\)', var_labels{var}, time(1), i, step_total),'Interpreter','latex');
end

% Animate over time
for i = 1:50:step_total
    for var = 1:7
        set(h(var), 'CData', squeeze(output_vars(var,:,:,1)));
        titles(var).String = ...
            sprintf('%s at $t=%.11f$ s \\(%d/%d\\)', var_labels{var}, time(i), i, step_total);
    end
    drawnow;
end

%% --- FUNCTIONS --- TODO: flux needs proper biasing in correctors
function [E, F] = flux(U, dx, dy, FD_method, R, cv, mu, k) 
    [~, nx, ny, ~] =  size(U);
    [E, F] = deal(zeros(4, nx,ny));
    [rho, u, v, T, p, ~, Et] = cons2prim(U, R, cv);
    if strcmpi(FD_method,'fwd') == 1 % predictor
        div_u = ddx_fwd(u,dx,1,true) + ddx_fwd(v,dy,2,true);
        tau_xy = mu.*(ddx_fwd(u,dy,2,true) + ddx_fwd(v,dx,1,true));
        tau_xx = 2*(mu.*ddx_fwd(u,dx,1,true) - (1/3)*div_u.*u);
        tau_yy = 2*(mu.*ddx_fwd(v,dy,2,true) - (1/3)*div_u.*v);
        qx = -k.*ddx_fwd(T,dx,1,true);
        qy = -k.*ddx_fwd(T,dy,2,true);
    elseif strcmpi(FD_method,'bwd') == 1 % predictor
        div_u = ddx_bwd(u,dx,1,true) + ddx_bwd(v,dy,2,true);
        tau_xy = mu.*(ddx_bwd(u,dy,2,true) + ddx_bwd(v,dx,1,true));
        tau_xx = 2*(mu.*ddx_bwd(u,dx,1,true) - (1/3).*div_u.*u);
        tau_yy = 2*(mu.*ddx_bwd(v,dy,2,true) - (1/3).*div_u.*v);
        qx = -k.*ddx_bwd(T,dx,1,true);
        qy = -k.*ddx_bwd(T,dy,2,true);
    elseif strcmpi(FD_method,'fwd_bias') == 1 % corrector
        div_u = ddx_fwd(u,dx,1,true) + ddx_fwd(v,dy,2,true);
        tau_xy = mu.*(ddx_central(u,dy,2,true) + ddx_bwd(v,dx,1,true));
        tau_xx = 2*(mu.*ddx_fwd(u,dx,1,true) - (1/3).*div_u.*u);
        tau_yy = 2*(mu.*ddx_fwd(v,dy,2,true) - (1/3).*div_u.*v);
        qx = -k.*ddx_fwd(T,dx,1,true);
        qy = -k.*ddx_fwd(T,dy,2,true);
    elseif strcmpi(FD_method,'bwd_bias') == 1 % corrector
        div_u = ddx_bwd(u,dx,1,true) + ddx_bwd(v,dy,2,true);
        tau_xy = mu.*(ddx_fwd(u,dy,2,true) + ddx_central(v,dx,1,true));
        tau_xx = 2*(mu.*ddx_bwd(u,dx,1,true) - (1/3).*div_u.*u);
        tau_yy = 2*(mu.*ddx_bwd(v,dy,2,true) - (1/3).*div_u.*v);
        qx = -k.*ddx_bwd(T,dx,1,true);
        qy = -k.*ddx_bwd(T,dy,2,true);
    else
        error("Invalid FD method given. Must be 'fwd', 'bwd', 'fwd_bias', 'bwd_bias'. ")
    end
    E(1,:,:) = rho.*u;
    E(2,:,:) = rho.*u.^2 + p - tau_xx;
    E(3,:,:) = rho.*u.*v - tau_xy;
    E(4,:,:) = (Et + p).*u - u.*tau_xx - v.*tau_xy + qx;
    F(1,:,:) = rho.*v;
    F(2,:,:) = rho.*u.*v - tau_xy;
    F(3,:,:) = rho.*v.^2 + p - tau_yy;
    F(4,:,:) = (Et + p).*v - v.*tau_yy - u.*tau_xy + qy;
end
function [p, u, v, T, rho, e, Et] = bc_enforcer(p, u, v, T, cv, R, uinf, p0, T0)
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


    rho = p./(R.*T);
    e = cv.*T;
    Et = rho .* (e + (u.^2 + v.^2)./2);
end