function fprime_y = ddy_fwd(f, dy, bc)
    fprime_y = zeros(size(f));
    [nx, ny] = size(f); 
    for i = 1:nx
        % fowards difference
        for j = 1:ny-1
            fprime_y(i,j) = (f(i,j+1) - f(i,j))/dy;
        end
        % boundary conditions
        if strcmpi(bc,'periodic') == 1
            fprime_y(i,ny) = (f(i,1) - f(i,ny))/dy;
        elseif strcmp(bc,'one-sided') == 1
            fprime_y(i,ny) = (f(i,ny) - f(i,ny-1))/dy;
        else
            error('No BC specified or specified incorrectly');
        end
    end
end
