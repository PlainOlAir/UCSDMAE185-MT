function fprime_x = ddx_fwd(f, dx,bc)
    fprime_x = zeros(size(f));
    [nx, ny] = size(f); 
    for j = 1:ny
        % forwards difference
        for i = 1:nx-1
            fprime_x(i,j) = (f(i+1,j) - f(i,j))/dx;
        end
        % boundary conditions
        if strcmpi(bc,'periodic') == 1
            fprime_x(nx,j) = (f(1,j) - f(nx,j))/dx;
        elseif strcmp(bc,'one-sided') == 1
            fprime_x(nx,j) = (f(nx,j) - f(nx-1,j))/dx;
        else
            error('No BC specified or specified incorrectly');
        end
    end
end
