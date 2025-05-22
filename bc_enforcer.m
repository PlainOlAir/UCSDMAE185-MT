function [rho, u, v, T, p, e, Et] = bc_enforcer(u, v, T, p, cv, R, uinf, pinf, Tinf)
    % Inlet (left boundary except corner)
    u(1, 2:end) = uinf;
    v(1, 2:end) = 0;
    p(1, 2:end) = pinf;
    T(1, 2:end) = Tinf;
    
    % Far field (top boundary)
    u(:, end) = uinf;
    v(:, end) = 0;
    p(:, end) = pinf;
    T(:, end) = Tinf;

    % Outlet (right boundary)
    u(end,:) = 2*u(end-1,:) - u(end-2,:);
    v(end,:) = 2*v(end-1,:) - v(end-2,:);
    p(end,:) = 2*p(end-1,:) - p(end-2,:);
    T(end,:) = 2*T(end-1,:) - T(end-2,:);

    % Bottom wall (adiabatic)
    u(2:end, 1) = 0;       % No-slip
    v(2:end, 1) = 0;       % No penetration
    % Adiabatic wall condition (dT/dy = 0)
    % T(2:end, 1) = T(2:end, 2);  % First-order approximation
    T(2:end,1) = (4*T(2:end,2) - T(2:end,3))/3; % 2nd order
    
    % Pressure at wall from momentum equation normal to wall
    p(2:end, 1) = p(2:end, 2);  % Simple approximation
    % Or more accurate: p(2:end,1) = (4*p(2:end,2) - p(2:end,3))/3;

    % Leading edge (corner point)
    u(1,1) = 0;
    v(1,1) = 0;
    p(1,1) = pinf;
    T(1,1) = Tinf;

    % Update dependent variables
    rho = p./(R.*T);
    e = cv.*T;
    Et = rho .* (e + (u.^2 + v.^2)./2);
end