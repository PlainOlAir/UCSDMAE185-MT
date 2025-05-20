clear
clc

addpath('lib\')

setup

dt = 2.35e-11; % Constant time step (s)
time = 0;

output_vars = cell(1,7);
convergence(:,:,1) = 1;
%% --- Main Loop ---
for step = 1:252
    % I/O, loop updates, delta_t_CFL, visualization
    a = sqrt(gamma*R*T);
    if ~isreal(T)
        disp(step)
        [row, col] = find(imag(T) ~= 0);
    end
    % rho, u, v, e, p, T, convergence
    % convergence_temp{step} = p_previous - p;
    convergence(:,:,step+1) = max(p_previous - p,[],'all');
    [rho, u, v, T, p, e, Et] = cons2prim(U, R, cv);
    % output_vars{1:6} = permute(cat(3, rho, u, v, e, p, T), [3 1 2]);
    new_vars = {rho, u, v, e, p, T, convergence};
    for k = 1:6
        if isempty(output_vars{k})
            output_vars{k} = new_vars{k};
        else
            output_vars{k} = cat(3, output_vars{k}, new_vars{k});
        end
    end
    output_vars{7} = convergence;
    % output_vars = {cat(3,output_vars{1}, rho), u, v, e, p, T};
    % output_vars{7} = convergence;
    p_previous = p;
    % compute delta_t CFL
    time(step+1) = time(step) + dt;
    
    %% Predictor    

    % update mu, k
    mu = sutherland(T);
    k = (cp/Pr)*mu;

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
    UBar = U - dt * (Edx + Fdy);

    % compute primitive var's from U_bar
    [rhoBar, uBar, vBar, TBar, pBar, eBar, EtBar] = cons2prim(UBar, R, cv); 

    % enforce BC's on p, u,v, T (update rho, e,...)
    [rhoBar, uBar, vBar, TBar, pBar, eBar, EtBar] = bc_enforcer(uBar, vBar, TBar, pBar, cv, R, uinf, pinf, Tinf);

    UBar = prim2cons(rhoBar, uBar, vBar, TBar, cv);

    % Corrector
    
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

end