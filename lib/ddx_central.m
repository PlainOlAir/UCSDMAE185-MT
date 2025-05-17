function firstdx = ddx_central(f, dx, periodicBC,  dir)
%% Description
% Calculates first derivative with second-order central difference 
% in the x direction

% INPUTS
% f = data matrix with the size of Nx x Ny
% dx = x-direction spacing of the data points
% periodicBC = true/false to set periodic boundaries (OPTIONAL)
% dir = direction (1 = x, 2 = y) (OPTIONAL)

% OUTPUTS
% firstdx = central FD first order derivative of input matrix f

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
firstdx = zeros(nx, ny); %allocate memory to output data


%% Calculate

% calculate boundary conditions depending on boundary type
if periodicBC %periodic boundary
    for j = 1:ny % forward FD
        firstdx(1,j) = (f(2,j) - f(end,j)) / (2 * dx);
        firstdx(end,j) = (f(1,j) - f(end-1,j)) / (2 * dx);
    end
else % non-periodic boudnary
    for j = 1:ny % forward & backwards FD at boundaries
        firstdx(1, j) = (-3 * f(1,j) + 4 * f(2,j) - f(3,j)) / (2 * dx);
        firstdx(end, j) = (3 * f(end,j) - 4 * f(end-1,j) + f(end-2,j)) / (2 * dx);
    end
end

for j = 1:ny %iterate through y points
    for i = 2:nx-1 %iterate through x points
        %central FD
        firstdx(i, j) = (f(i+1,j) - f(i-1,j)) / (2 * dx);
    end
end

if dir == 2 % untranspose data matrix if transposed earlier
    firstdx = firstdx';
end

end