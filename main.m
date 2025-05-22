clear
clc

addpath('lib\')

setup

time = 0;

output_vars = cell(1,7);
for i = 1:6
    output_vars{i} = zeros(nx, ny, step_total+1);
    output_vars{i}(:, :, 1) = new_vars{i};
end

convergence = zeros(step_total+1, 4);

%% --- Main Loop ---
for step = 1:step_total
    % I/O, loop updates, delta_t_CFL, visualization

    %%%%%% Predictor %%%%%%

    % update mu, k
    mu = sutherland(T);
    k = (cp/Pr)*mu;
    
    % compute delta_t CFL
    a = sqrt(gamma .* R .* T);
    vprime = max(4 / 3 .* mu .* (k .* mu ./ Pr) ./ rho, [], 'all');
    dtCFL = (abs(u)./dx + abs(v)./dy + a .* sqrt(1/(dx^2) + 1/(dy^2)) + 2 .* vprime .* (1/(dx^2) + 1/(dy^2))).^(-1);
    dt = 2.35e-11; % Constant time step (s)
    time(step+1) = time(step) + dt;
    %dt = permute(repmat(dtCFL,[1 1 4]), [3 1 2]);
    
    % compute derivatives, update E, F
    if mod(step, 2) == 0
        [E, F] = flux('backward', U, dx, dy, mu, k, R, cv);
        Edx = dxn_3d(@ddx_fwd, E, dx);
        Fdy = dxn_3d(@ddx_fwd, F, dy, 2);
    else
        [E, F] = flux('forward', U, dx, dy, mu, k, R, cv);
        Edx = dxn_3d(@ddx_bwd, E, dx);
        Fdy = dxn_3d(@ddx_bwd, F, dy, 2);
    end

    % compute U_bar using U, E, F
    UBar = U - dt .* (Edx + Fdy);

    % compute primitive var's from U_bar
    [rhoBar, uBar, vBar, TBar, pBar, eBar, EtBar] = cons2prim(UBar, R, cv);

    % enforce BC's on p, u,v, T (update rho, e,...)
    [rhoBar, uBar, vBar, TBar, pBar, eBar, EtBar] = bc_enforcer(uBar, vBar, TBar, pBar, cv, R, uinf, pinf, Tinf);

    UBar = prim2cons(rhoBar, uBar, vBar, TBar, cv);

    %%%%% Corrector %%%%%

    % update mu, k
    muBar = sutherland(TBar);
    kBar = (cp/Pr)*muBar;

    % compute derivatives, update E, F
    if mod(step, 2) == 0
        [EBar, FBar] = flux('forward', UBar, dx, dy, muBar, kBar, R, cv);
        EBardx = dxn_3d(@ddx_bwd, EBar, dx);
        FBardy = dxn_3d(@ddx_bwd, FBar, dy, 2);
    else
        [EBar, FBar] = flux('backward', UBar, dx, dy, muBar, kBar, R, cv);
        EBardx = dxn_3d(@ddx_fwd, EBar, dx);
        FBardy = dxn_3d(@ddx_fwd, FBar, dy, 2);
    end

    % compute U from primitive vars
    U = 0.5 .* (U + UBar - dt.*(EBardx + FBardy));

    % compute primitive vars from U
    [rho, u, v, T, p, e, Et] = cons2prim(U, R, cv);

    % enforce BC's on p, u, v, T (update rho, e,...)
    [rho, u, v, T, p, e, Et] = bc_enforcer(u, v, T, p, cv, R, uinf, pinf, Tinf);

    % compute U from primitive vars
    U = prim2cons(rho,u,v,T,cv);

    %%%%% Data Variables %%%%%%

    [rho, u, v, T, p, e, Et] = cons2prim(U, R, cv);

    new_vars = {rho, u, v, e, p, T};

    for i = 1:4
        % convergence(step+1, i) = max(abs(U(i, :, :) - U_prev(i, :, :)) ./ (0.5 * (abs(U(i, :, :)) + abs(U_prev(i, :, :))) + eps), [], 'all');
        % convergence(step+1, i) = sqrt(mean((U(i,:,:) - U_prev(i,:,:)).^2, 'all'));
        convergence(step+1, i) = sqrt(mean((U(i,:,:) - U_prev(i,:,:)).^2, 'all')) / max(abs(U(i,:,:)), [], 'all');
    end

    for k = 1:6
        output_vars{k}(:, :, step) = new_vars{k};
    end
    U_prev = U;
end
output_vars{7} = convergence;