clear

addpath('lib\')

setup

dt = 2.35e-11; % Constant time step (s)

time = 0;

%% --- Main Loop ---
for step = 1:step_total
    %% I/O, loop updates, delta_t_CFL, visualization
    a = sqrt(gamma*R*T);
    % compute delta_t CFL
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
    [E, F] = flux(U,dx,dy,FD_method_pred);
    % compute U_bar using U, E, F
    U_bar = U(:,:,step) - dt.*(E(:,:,step) + F(:,:,step));
    % compute primitive var's from U_bar
    [rho_bar, u_bar, v_bar, T_bar, p_bar, e_bar, Et_bar] = cons2prim(U_bar, R, cv); 
    % enforce BC's on p, u,v, T (update rho, e,...)


    %% Corrector
    [E, F] = flux(U_pred,dx,dy,FD_method_corr);

    %% Compute U from primitive vars

end
%% --- Animator ---
%% Draw figure
figure;
tile = tiledlayout(1,1);
axesArray = [nexttile(tile,1)];
h = [pcolor(axesArray(1),xx, yy, squeeze(U(:,:,1)))];
for ax = axesArray
    shading(ax, 'interp');
    axis(ax, 'equal', 'tight');
    xlabel(ax, 'x');
    ylabel(ax, 'y');
    colorbar(ax);
    % clim(ax,[0 2])
end
titles = [title(axesArray(1),'', 'Interpreter', 'latex')];
for i = 1:10:step_total
    set(h(1), 'CData', squeeze(U(:,:,i)));
    titles(1).String = sprintf('$\\ U \\ at \\ t=%.3f s \\ (%d/%d)$',t(1,i), i, step_total);  % Update title
    drawnow;
end

%% --- FUNCTIONS ---
function [E, F] = flux(U, dx, dy, FD_method) 
% need gradient and bias eq's
    [rho, u, v, T, p, ~, Et] = cons2prim(U, R, cv);
    if strcmpi(FD_method,'fwd') == 1 % predictor
        tau_xy = mu(ddx_fwd(u,dy,2,true) + ddx_fwd(v,dx,1,true));
        tau_xx = 2.*(mu.*ddx_fwd(u,dx,1,true) - (1/3).*gradient.*u);
        tau_yy = 2.*(mu.*ddx_fwd(v,dy,2,true) - (1/3).*gradient.*v);
        qx = -k.*ddx_fwd(T,dx,1,true);
        qy = -k.*ddx_fwd(T,dy,2,true);
    elseif strcmpi(FD_method,'bwd') == 1 % predictor
        tau_xy = mu(ddx_bwd(u,dy,2,true) + ddx_bwd(v,dx,1,true));
        tau_xx = 2.*(mu.*ddx_bwd(u,dx,1,true) - (1/3).*gradient.*u);
        tau_yy = 2.*(mu.*ddx_bwd(v,dy,2,true) - (1/3).*gradient.*v);
        qx = -k.*ddx_bwd(T,dx,1,true);
        qy = -k.*ddx_bwd(T,dy,2,true);
    elseif strcmpi(FD_method,'fwd_bias') == 1 % corrector
        tau_xy = mu(ddx_central(u,dy,2,true) + ddx_bwd(v,dx,1,true));
        tau_xx = 2.*(mu.*ddx_fwd(u,dx,1,true) - (1/3).*gradient.*u);
        tau_yy = 2.*(mu.*ddx_fwd(v,dy,2,true) - (1/3).*gradient.*v);
        qx = -k.*ddx_fwd(T,dx,1,true);
        qy = -k.*ddx_fwd(T,dy,2,true);
    elseif strcmpi(FD_method,'bwd_bias') == 1 % corrector
        tau_xy = mu(ddx_fwd(u,dy,2,true) + ddx_central(v,dx,1,true));
        tau_xx = 2.*(mu.*ddx_bwd(u,dx,1,true) - (1/3).*gradient.*u);
        tau_yy = 2.*(mu.*ddx_bwd(v,dy,2,true) - (1/3).*gradient.*v);
        qx = -k.*ddx_bwd(T,dx,1,true);
        qy = -k.*ddx_bwd(T,dy,2,true);
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