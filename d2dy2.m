function f2prime_y = d2dy2(f, dy,bc)
    f2prime_y = zeros(size(f));
    [nx, ny] = size(f);
    for i = 1:nx
        % central difference
        for j = 2:ny-1
            f2prime_y(i,j) = (f(i,j+1) -2*f(i,j) + f(i,j-1))/dy^2;
        end
        % boundary conditions
        if strcmpi(bc,'periodic') == 1
            f2prime_y(i,1) = (f(i,2) -2*f(i,1) +f(i,ny))/dy^2;
            f2prime_y(i,ny) = (f(i,1) -2*f(i,ny) +f(i,ny-1))/dy^2;
        elseif strcmp(bc,'one-sided') == 1
            f2prime_y(i,1) = (2*f(i,1) -5*f(i,2) +4*f(i,3) -f(i,4))/dy^2;
            f2prime_y(i,ny) = (2*f(i,ny) -5*f(i,ny-1) +4*f(i,ny-2) -f(i,ny-3))/dy^2;
        else
            error('No BC specified or specified incorrectly');
        end
    end
end