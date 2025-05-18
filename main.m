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

end

function F = yflux(U, p, T, cp, Pr, dx, dy)
end

function pred = predictor()
    [rho, u, v, T, p, e, Et] = cons2prim(U, R, cv);
    % update mu, k
    mu = sutherland(T, mu0, T0, S1);
    % k = ???;
    % compute derivatives
    tau_xy = mu(ddx_fwd(u,dy,2,true) + ddx_fwd(v,dx,1,true));
    tau_xx = 2.*(mu.*ddx_fwd(u,dx,1,true) - (1/3).*gradient.*u);
    tau_yy = 2.*(mu.*ddx_fwd(v,dy,2,true) - (1/3).*gradient.*v);
    qx = -k.*ddx_fwd(T,dx,1,true);
    qy = -k.*ddx_fwd(T,dy,2,true);
            % fwd diff1 for dE/dx
    E(1,1) = rho.*u;
    E(2,1) = rho.*u.^2 + p - tau_xx;
    E(3,1) = rho.*u.*v - tau_xy;
    E(4,1) = (Et + p).*u - u.*tau_xx - v.*tau_xy + qx;
        % fwd diff1 for dF/dy
    F(1,1) = rho.*v;
    F(2,1) = rho.*u.*v - tau_xy;
    F(3,1) = rho.*v.^2 + p - tau_yy;
    F(4,1) = (Et + p).*v - v.*tau_yy - u.*tau_xy + qy;


    dE_dx = d2dx2(E,dx,1,true);
    dE_dy = d2dx2(E,dy,2,true);
    E_pred = 
    pred = U - dt*(ddx_fwd(u(:,:,step),dx,'periodic') + ddy_fwd(v(:,:,step),dy,'periodic'));
    % update E

    % 
end

function corr = corrector()
    % bwd diff1 for dE/dx
    % bwd diff1 for dE/dy
    % Fbar bias: tau_xy = mu(fwd_diff1y + cent_diff2x
end