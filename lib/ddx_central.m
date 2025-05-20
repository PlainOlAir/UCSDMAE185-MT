function firstdx = ddx_central(f, dx, dir, periodicBC)
%% Description
% Calculates first derivative with second-order central difference
% in the x direction

% INPUTS
% f = data matrix with the size of Nx x Ny
% dx = x-direction spacing of the data points
% dir = direction (1 = x, 2 = y) (OPTIONAL)
% periodicBC = true/false to set periodic boundaries (OPTIONAL)

% OUTPUTS
% firstdx = central FD first order derivative of input matrix f

%% Setup
% process input variables
switch nargin
    case 2
        dir = 1;
        periodicBC = false;
    case 3
        periodicBC = false;
end

if dir == 2 % if going in y direction, transpose data matrix
    f = f';
end

[nx, ny] = size(f); %get number of x and y points
firstdx = zeros(nx, ny); %allocate memory to output data


%% Calculate

% calculate boundary conditions depending on boundary type
if periodicBC %periodic boundary
    firstdx(1,:) = (f(2,:) - f(end,:)) / (2 * dx);
    firstdx(end,:) = (f(1,:) - f(end-1,:)) / (2 * dx);
else % non-periodic boudnary
    firstdx(1, :) = (-3 * f(1,:) + 4 * f(2,:) - f(3,:)) / (2 * dx);
    firstdx(end, :) = (3 * f(end,:) - 4 * f(end-1,:) + f(end-2,:)) / (2 * dx);
end

for i = 2:nx-1 %iterate through x points
    %central FD
    firstdx(i, :) = (f(i+1,:) - f(i-1,:)) / (2 * dx);
end

if dir == 2 % untranspose data matrix if transposed earlier
    firstdx = firstdx';
end

end