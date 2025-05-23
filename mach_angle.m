%%%%%%%%%%%%%%%%%%
%%% Mach Angle %%%
%%%%%%%%%%%%%%%%%%

figure; 

theta = asin(1/M);
y_end = L * tan(theta);

h = pcolor(xx, yy, output_vars{3}(:, :, end));
hold on
h2 = plot([0 L], [0 y_end], 'LineWidth', 2, 'Color', 'red');
hold off

shading interp
axis equal tight
xlabel('$x ~ (m)$', 'Interpreter', 'latex');
ylabel('$y ~ (m)$', 'Interpreter', 'latex');
colorbar;
title(sprintf('M = %d', M), 'Interpreter', 'latex');
legend({'', 'Theoretical Mach Angle'})