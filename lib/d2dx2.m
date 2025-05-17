function scnddx = d2dx2(f, dx, periodicBC, dir)
%% Description
% Calculates second derivative with second-order central difference
% in the x direction

% INPUTS
% f = data matrix with the size of Nx x Ny
% dx = x-direction spacing of the data points
% periodicBC = true/false to set periodic boundaries (OPTIONAL)
% dir = direction (1 = x, 2 = y) (OPTIONAL)

% OUTPUTS
% scnddy = central FD second order derivative of input matrix f

%% Setup
% process input variables
switch nargin
    case 2
        dir = 1;
        periodicBC = false;
    case 3
        dir = 1;
end

if dir == 2 % if going in y direction, transpose data matrix
    f = f';
end

[nx, ny] = size(f); %get number of x and y points
scnddx = zeros(nx, ny); %allocate memory to output data
%% Calculate

% calculate boundary conditions depending on boundary type
if periodicBC % periodic boundary
    for j = 1:ny
        scnddx(1,j) = (f(2,j) - 2 * f(1,j) + f(end,j)) / (dx^2);
        scnddx(end,j) = (f(1,j) - 2 * f(end,j) + f(end-1,j)) / (dx^2);
    end
else % non-periodic boundary
    for j = 1:ny
        scnddx(1,j) = (2 * f(1,j) - 5 * f(2,j) + 4 * f(3,j) - f(4,j)) / (dx^2); % forward FD
        scnddx(end,j) = (2 * f(end,j) - 5 * f(end-1,j) + 4 * f(end-2,j) - f(end-3,j)) / (dx^2); %backwards FD
    end
end

% calculate central derivatives
for j = 1:ny %iterate through y points
    for i = 2:nx-1 %iterate through x points
        %central FD
        scnddx(i,j) = (f(i+1,j) - 2 * f(i,j) + f(i-1,j)) / (dx^2);
    end
end

if dir == 2 % untranspose data matrix if transposed earlier
    scnddx = scnddx';
end

end