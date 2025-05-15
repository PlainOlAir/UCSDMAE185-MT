function fprime_x = ddx_bwd(f, dx,bc)
    fprime_x = zeros(size(f));
    [nx, ny] = size(f); 
    for j = 1:ny
        % backwards difference
        for i = 2:nx
            fprime_x(i,j) = (f(i,j) - f(i-1,j))/dx;
        end
        % boundary conditions
        if strcmpi(bc,'periodic') == 1
            fprime_x(1,j) = (f(1,j) - f(nx,j))/dx;
        elseif strcmp(bc,'one-sided') == 1
            fprime_x(1,j) = (f(2,j) - f(1,j))/dx;
        else
            error('No BC specified or specified incorrectly');
        end
    end
end
