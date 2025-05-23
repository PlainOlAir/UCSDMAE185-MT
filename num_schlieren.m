%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Numerical Schlieren %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Setup

% Calculate refractive index from density gradients
gradp = sqrt(ddx_central(output_vars{5}(:,:,step_total), dx, 1) .^ 2 + ddx_central(output_vars{5}(:,:,step_total), dy, 2) .^ 2);

% Fudge factors
beta = 0.8;
kappa = 10;

% Apply fudge factors (makes a Schlieren-like image)
S = beta .* exp(-kappa .* (gradp ./ (max(gradp))));

%% Plotting
figure;
pcolor(xx,yy,S)
shading('interp')
colormap(gray)
colorbar
clim([0, 1])
xlabel('$x ~ (m)$', 'Interpreter','latex')
ylabel('$y ~ (m)$', 'Interpreter','latex')
title('Numerical Schlieren','Interpreter','latex')
axis equal tight


