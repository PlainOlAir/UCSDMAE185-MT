clear
close all
clc

addpath('lib\')

setup

dt = 2.35e-11; % Constant time step (s)
time = 0;

%% --- Main Loop ---
for step = 1:step_total
    % I/O, loop updates, delta_t_CFL, visualization
    a = sqrt(gamma*R*T);

    % rho, u, v, e, p, T, convergence
    convergence = p_previous - p;
    [rho, u, v, T, p, e, Et] = cons2prim(U(:,:,:), R, cv);
    output_vars(:, :, :, step) = permute(cat(3, rho, u, v, e, p, T, convergence), [3 1 2]);
    p_previous = p;
    % compute delta_t CFL
    time(step+1) = time(step) + dt;

    %% Predictor

    % update mu, k
    mu = sutherland(T);
    k = (cp/Pr)*mu;

    % compute derivatives, update E, F
    if mod(step, 2) == 1
        [E, F] = flux('backward', U(:,:,:), dx, dy, mu, k, R, cv);
        Edx = dxn_3d(@ddx_fwd, E, dx);
        Fdy = dxn_3d(@ddx_fwd, F, dy, 2);
    else
        [E, F] = flux('forward', U(:,:,:), dx, dy, mu, k, R, cv);
        Edx = dxn_3d(@ddx_bwd, E, dx);
        Fdy = dxn_3d(@ddx_bwd, F, dy, 2);
    end

    % compute U_bar using U, E, F
    UBar = U(:,:,:) - dt * (Edx + Fdy);

    % compute primitive var's from U_bar
    [rhoBar, uBar, vBar, TBar, pBar, eBar, EtBar] = cons2prim(UBar, R, cv); 

    % enforce BC's on p, u,v, T (update rho, e,...)
    [rhoBar, uBar, vBar, TBar, pBar, eBar, EtBar] = bc_enforcer(uBar, vBar, TBar, pBar, cv, R, uinf, p0, T0);

    UBar = prim2cons(rhoBar, uBar, vBar, TBar, cv);

    % Corrector
    
    % update mu, k
    muBar = sutherland(TBar);
    kBar = (cp/Pr)*muBar;

    % compute derivatives, update E, F
    if mod(step, 2) == 1
        [EBar, FBar] = flux('forward', U(:,:,:), dx, dy, mu, k, R, cv);
        EBardx = dxn_3d(@ddx_bwd, EBar, dx);
        FBardy = dxn_3d(@ddx_bwd, FBar, dy, 2);
    else
        [EBar, FBar] = flux('backward', U(:,:,:), dx, dy, mu, k, R, cv);
        EBardx = dxn_3d(@ddx_fwd, EBar, dx);
        FBardy = dxn_3d(@ddx_fwd, FBar, dy, 2);
    end

    % compute U from primitive vars
    U(:,:,:) = 0.5.*(U(:,:,:) + UBar - dt.*(EBardx + FBardy));

    % compute primitive vars from U
    [rho, u, v, T, p, e, Et] = cons2prim(U(:,:,:), R, cv);

    % enforce BC's on p, u, v, T (update rho, e,...)
    [rho, u, v, T, p, e, Et] = bc_enforcer(u, v, T, p, cv, R, uinf, p0, T0);

    % compute U from primitive vars
    U(:,:,:) = prim2cons(rho,u,v,T,cv);

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
    colorbar(axesArray(var));
    titles(var) = title(axesArray(var), ...
        sprintf('%s at $t=%.11f$ s \\(%d/%d\\)', var_labels{var}, time(1), i, step_total),'Interpreter','latex');
end
% Animate over time
for i = 2:50:step_total
    for var = 1:7
        h(var).CData = squeeze(output_vars(var, :, :, i));
        titles(var).String = ...sprintf('%s at $t=%.11f$ s \\(%d/%d\\)', var_labels{var}, time(i), i, step_total);
    end
    drawnow;
end