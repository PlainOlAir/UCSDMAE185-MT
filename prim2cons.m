function [U] = prim2cons(rho,u,v,T,cv)
[nx,ny]     = size(rho);          
U           = zeros(4,nx,ny);

U(1,:,:)	= rho;
U(2,:,:)    = rho.*u;
U(3,:,:)    = rho.*v;
U(4,:,:)    = rho.*(cv*T + (u.^2+v.^2)/2);
end