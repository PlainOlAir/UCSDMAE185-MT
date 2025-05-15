function fprime_y = ddy_bwd(f, dy, bc)
    fprime_y = zeros(size(f));
    [nx, ny] = size(f); 
    for i = 1:nx
        % backwards difference
        for j = 2:ny
            fprime_y(i,j) = (f(i,j) - f(i,j-1))/dy;
        end
        % boundary conditions
        if strcmpi(bc,'periodic') == 1
            fprime_y(i,1) = (f(i,1) - f(i,ny))/dy;
        elseif strcmp(bc,'one-sided') == 1
            fprime_y(i,j) = (f(i,2) - f(i,1))/dy;
        else
            error('No BC specified or specified incorrectly');
        end
    end
end
