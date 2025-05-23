function [E, F] = flux(FD_method, U, dx, dy, mu, k, R, cv) 
    [~, nx, ny] = size(U);
    [E, F] = deal(zeros(4, nx, ny));
    [rho, u, v, T, p, ~, Et] = cons2prim(U, R, cv);

    % Compute derivatives based on FD method
    if strcmpi(FD_method,'forward')
        div_u = ddx_fwd(u, dx, 1) + ddx_fwd(v, dy, 2);
        tau_xyE = mu .* (ddx_central(u, dy, 2) + ddx_fwd(v, dx, 1));
        tau_xyF = mu .* (ddx_fwd(u, dy, 2) + ddx_central(v, dx, 1));
        tau_xx = 2 * mu .* (ddx_fwd(u, dx, 1) - (1/3) .* div_u);
        tau_yy = 2 * mu .* (ddx_fwd(v, dy, 2) - (1/3) .* div_u);
        qx = -k .* ddx_fwd(T, dx, 1);
        qy = -k .* ddx_fwd(T, dy, 2);
    else
        div_u = ddx_bwd(u, dx, 1) + ddx_bwd(v, dy, 2);
        tau_xyE = mu .* (ddx_central(u, dy, 2) + ddx_bwd(v, dx, 1));
        tau_xyF = mu .* (ddx_bwd(u, dy, 2) + ddx_central(v, dx, 1));
        tau_xx = 2 .* mu .* (ddx_bwd(u, dx, 1) - (1/3) .* div_u);
        tau_yy = 2 .* mu .* (ddx_bwd(v, dy, 2) - (1/3) .* div_u);
        qx = -k .* ddx_bwd(T, dx, 1);
        qy = -k .* ddx_bwd(T, dy, 2);
    end
    
    % Compute regular fluxes
    E(1,:,:) = rho .* u;
    E(2,:,:) = rho .* u.^2 + p - tau_xx;
    E(3,:,:) = rho .* u .* v - tau_xyE;
    E(4,:,:) = (Et + p) .* u - u .* tau_xx - v .* tau_xyE + qx;

    F(1,:,:) = rho .* v;
    F(2,:,:) = rho .* u .* v - tau_xyF;
    F(3,:,:) = rho .* v.^2 + p - tau_yy;
    F(4,:,:) = (Et + p) .* v - v .* tau_yy - u .* tau_xyF + qy;
end