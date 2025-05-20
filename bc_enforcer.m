function [rho, u, v, T, p, e, Et] = bc_enforcer(u, v, T, p, cv, R, uinf, p0, T0)
    % Leading edge
    u(1,1) = 0;
    v(1,1) = 0;
    p(1,1) = p0;
    T(1,1) = T0;

    % Bottom wall
    u(:, 1) = 0;
    v(:, 1) = 0;
    T(:, 1) = T0;
    p(:, 1) = 2*(p(:,3)) - p(:,2);
    
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
    
    % Outlet
    u(end,:) = 2*(u(end-2,:)) - u(end-1,:);
    v(end,:) = 2*(v(end-2,:)) - v(end-1,:);
    p(end,:) = 2*(p(end-2,:)) - p(end-1,:);
    T(end,:) = 2*(T(end-2,:)) - T(end-1,:);

    rho = p./(R.*T);
    e = cv.*T;
    Et = rho .* (e + (u.^2 + v.^2)./2);
end