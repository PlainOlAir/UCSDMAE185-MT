function [rho, u, v, T, p, e, Et] = bc_enforcer(u, v, T, p, cv, R, uinf, pinf, Tinf)
    % Outlet
    u(end,:) = 2*(u(end-2,:)) - u(end-1,:);
    v(end,:) = 2*(v(end-2,:)) - v(end-1,:);
    p(end,:) = 2*(p(end-2,:)) - p(end-1,:);
    T(end,:) = 2*(T(end-2,:)) - T(end-1,:);

    % Bottom wall
    u(:, 1) = 0;
    v(:, 1) = 0;
    T(:, 1) = Tinf;
    p(:, 1) = 2*(p(:,3)) - p(:,2);

    % Leading edge
    u(1,1) = 0;
    v(1,1) = 0;
    p(1,1) = pinf;
    T(1,1) = Tinf;
    
    % Inlet
    u(1, :) = uinf;
    v(1, :) = 0;
    p(1, :) = pinf;
    T(1, :) = Tinf;
    
    % Far field
    u(:, end) = uinf;
    v(:, end) = 0;
    p(:, end) = pinf;
    T(:, end) = Tinf;

    rho = p./(R.*T);
    e = cv.*T;
    Et = rho .* (e + (u.^2 + v.^2)./2);
end