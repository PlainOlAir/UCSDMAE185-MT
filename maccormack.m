%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MacCormack Solver %%%
%%%%%%%%%%%%%%%%%%%%%%%%%

%% Description %%
% Applies the MacCormack method to solve the 2-D Navier-Stokes compressible
% flows
% Compact vector form:
% dU/dt + dE/dx + dF/dy = 0

for step = 1:step_total
    %% MacCormack Loop
    
    % Calculate mu & k with new U
    mu = sutherland(T);
    k = (cp/Pr)*mu;
    
    % Compute dt with CFL criterion
    a = sqrt(gamma .* R .* T);
    vprime = mu ./ rho + k ./ (rho * cp);
    dt = min((abs(u)./dx + abs(v)./dy + a .* sqrt(1/(dx^2) + 1/(dy^2)) + 2 .* vprime .* (1/(dx^2) + 1/(dy^2))).^(-1), [], 'all');

    % Add time step
    time(step+1) = time(step) + dt; 

    % Build as a 4 X nx X ny array for multiplication
    dt = permute(repmat(dt,[1 1 4]), [3 1 2]); 

    %% Predictor Step %%
    
    % Compute E & F and apply derivatives
    % Alternate between backward and forward FD stencils per time step
    if mod(step, 2) == 0

        % Flux calculation
        [E, F] = flux('backward', U, dx, dy, mu, k, R, cv);

        % Apply derivatives
        Edx = dxn_3d(@ddx_fwd, E, dx);
        Fdy = dxn_3d(@ddx_fwd, F, dy, 2);

    else

        % Flux calculation
        [E, F] = flux('forward', U, dx, dy, mu, k, R, cv);

        % Apply derivatives
        Edx = dxn_3d(@ddx_bwd, E, dx);
        Fdy = dxn_3d(@ddx_bwd, F, dy, 2);
        
    end

    % Compute UBar
    UBar = U - dt .* (Edx + Fdy);

    % Compute primitive variables from UBar
    [rhoBar, uBar, vBar, TBar, pBar, eBar, EtBar] = cons2prim(UBar, R, cv, nx, ny);

    % Enforce BCs
    [rhoBar, uBar, vBar, TBar, pBar, eBar, EtBar] = bc_enforcer(uBar, vBar, TBar, pBar, cv, R, uinf, pinf, Tinf);
    
    % Rebuild UBar with BCs
    UBar = prim2cons(rhoBar, uBar, vBar, TBar, cv);

    % Calculate muBar and kBar with predictor step values
    muBar = sutherland(TBar);
    kBar = (cp/Pr)*muBar;

    %% Corrector Step %%

    % Compute EBar & FBar and apply derivatives
    % Alternate between backward and forward FD stencils per time step
    % OPPOSITE OF PREDICTOR DIRECTION
    if mod(step, 2) == 0

        % Flux calculation
        [EBar, FBar] = flux('forward', UBar, dx, dy, muBar, kBar, R, cv);

        % Apply derivatives
        EBardx = dxn_3d(@ddx_bwd, EBar, dx);
        FBardy = dxn_3d(@ddx_bwd, FBar, dy, 2);

    else
        
        % Flux calculation
        [EBar, FBar] = flux('backward', UBar, dx, dy, muBar, kBar, R, cv);

        % Apply derivatives
        EBardx = dxn_3d(@ddx_fwd, EBar, dx);
        FBardy = dxn_3d(@ddx_fwd, FBar, dy, 2);

    end

    % Solve for U, applying both predictor and corrector steps
    U = 0.5 .* (U + UBar - dt.*(EBardx + FBardy));

    % Compute primitive variables from new U
    [rho, u, v, T, p, e, Et] = cons2prim(U, R, cv, nx, ny);

    % Enforce BCs
    [rho, u, v, T, p, e, Et] = bc_enforcer(u, v, T, p, cv, R, uinf, pinf, Tinf);

    % Rebuild U with BCs
    U = prim2cons(rho,u,v,T,cv);

    %% Data Processing %%
    
    % Build array of values to be stored
    new_vars = {rho, u, v, e, p, T};

    % Store current step values
    for k = 1:6
        output_vars{k}(:, :, step+1) = new_vars{k};
    end
    
    % Calculate normalized residuals
    output_vars{7}(step+1, :) = sqrt(mean((U - U_prev).^2, [2 3])) ./ max(abs(U), [], [2 3]);
    U_prev = U;
end