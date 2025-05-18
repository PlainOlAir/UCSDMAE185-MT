clear

addpath('lib\')

setup

dt = 2.35e-11; % Constant time step (s)

time = 0;

for step = 1:1500
    % primitive var's and U update
    U = prim2cons(rho(:, :, step), u(:, :, step), v(:, :, step), T(:, :, step), cv);
    % I/O and visualization
    % compute delta_t CFL
    %% Predictor
    if mod(step, 2) == 1
        
    else
    end

    if mod(step, 2) == 1

    else
    end
end

function E = xflux(U, p, T, cp, Pr, dx, dy)
    % tau_xy = mu(du/dy + dv/dx))
    % tau_xx = 2(mu*du/dx - 1/3*gradient.u)
    % tau_yy = 2(mu*dv/dy - 1/3*gradient.v)
    % qx = -k*dT/dx
    % qy = -k*dT/dy
    [rho, u, v, T, p, e, Et] = cons2prim(U, R, cv);
    E(1,1) = U(2,1);
    E(2,1) = (U(2,1)^2)/rho + p - tau_xx;
    E(3,1) = U(2,1)*U(3,1)/rho - tau_xy;
    E(4,1) = (U(4,1) + p)*U(2,1)/rho - U(1,1)
end

function F = yflux(U, p, T, cp, Pr, dx, dy)
end

function pred = predictor()
    % update mu, k
    mu = sutherland(T, mu0, T0, S1);
    % k = ???;
    % compute derivatives
    dE_dx = d2dx2(E,dx,1,true);
    dE_dy = d2dx2(E,dy,2,true);
    E_pred = 
    pred = U - dt*(ddx_fwd(u(:,:,step),dx,'periodic') + ddy_fwd(v(:,:,step),dy,'periodic'));
    % update E
    % 
end