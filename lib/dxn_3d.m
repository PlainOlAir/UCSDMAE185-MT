function dx3 = dxn_3d(func, f, dx, dir, periodicBC)

switch nargin
    case 3
        dir = 1;
        periodicBC = false;
    case 4
        periodicBC = false;
end

[nvars, nx, ny] = size(f);
dx3 = zeros(nvars, nx, ny);

for i = 1:nvars
    dx3(i, :, :) = func(squeeze(f(i, :, :)), dx, dir, periodicBC);
end

end