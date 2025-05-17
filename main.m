clear

addpath('lib\')

setup

dt = 2.35e-11; % Constant time step (s)

time = 0;

for step = 1:1500
    U = prim2cons(rho(:, :, step), u(:, :, step), v(:, :, step), T(:, :, step), cv);
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