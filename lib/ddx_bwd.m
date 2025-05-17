function firstdx = ddx_bwd(f, dx, dir, periodicBC)
%% Description
% Calculates first derivative with first-order backwards difference 
% in the x direction

% INPUTS
% f = data matrix with the size of Nx x Ny
% dx = x-direction spacing of the data points
% dir = direction (1 = x, 2 = y) (OPTIONAL)
% periodicBC = true/false to set periodic boundaries (OPTIONAL)

% OUTPUTS
% firstdx = backwards FD first derivative of input matrix f

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
    for j = 1:ny % backwards FD
        firstdx(1, j) = (f(1, j) - f(end,j)) / dx;
    end
else % non-periodic boudnary
    for j = 1:ny % forward FD
        firstdx(1, j) = (f(2,j) - f(1,j)) / dx;
    end
end

for j = 1:ny %iterate through y points
    for i = 2:nx %iterate through x points
        %backwards FD
        firstdx(i, j) = (f(i,j) - f(i-1,j)) / dx;
    end
end

if dir == 2 % untranspose data matrix if transposed earlier
    firstdx = firstdx';
end

end
