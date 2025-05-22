figure;

beta = 0.8;
kappa = 10;
gradp = sqrt(ddx_central(output_vars{5}(:,:,step_total),dx,1).^2 + ddx_central(output_vars{5}(:,:,step_total),dy,2).^2);
S = beta*exp(-kappa*(gradp./(max(gradp))));
pcolor(xx,yy,S)
shading('interp')
colormap(gray)
colorbar
clim([0, 1])
xlabel('$x$', 'Interpreter','latex')
ylabel('$y$', 'Interpreter','latex')
title('Numerical Schlieren','Interpreter','latex')
axis equal tight


