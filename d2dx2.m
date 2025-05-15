function f2prime_x = d2dx2(f, dx,bc)
    f2prime_x = zeros(size(f));
    [nx, ny] = size(f);
    for j = 1:ny
        % central difference
        for i = 2:nx-1
            f2prime_x(i,j) = (f(i+1,j) -2*f(i,j) + f(i-1,j))/dx^2;
        end
        % boundary conditions
        if strcmpi(bc,'periodic') == 1
            f2prime_x(1,j) = (f(2,j) -2*f(1,j) +f(nx,j))/dx^2;
            f2prime_x(nx,j) = (f(1,j) -2*f(nx,j) +f(nx-1,j))/dx^2;
        elseif strcmp(bc,'one-sided') == 1
            f2prime_x(1,j) = (2*f(1,j) -5*f(2,j) +4*f(3,j) -f(4,j))/dx^2;
            f2prime_x(nx,j) = (2*f(nx,j) -5*f(nx-1,j) +4*f(nx-2,j) -f(nx-3,j))/dx^2;
        else
            error('No BC specified or specified incorrectly');
        end
    end
end