length = sqrt((1e-5)^2 + (H)^2); 
theta = asind(1/M)
[x1,y1] = deal([0, 0]);
x_end = length * cosd(theta);
y_end = length * sind(theta);
tile = tiledlayout(1, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
% initialize handles
axesArray = gobjects(1, 1);
h = gobjects(1, 1);
titles = gobjects(1, 1);
for var = 1:1
    axesArray(var) = nexttile(tile, var);
    output_frame = squeeze(output_vars{var}(:,:,step_total));
    hold on
    h(var) = pcolor(axesArray(var), xx, yy, output_frame);
    h2(var) = line(axesArray(var), [0 x_end], [0 y_end]);
    hold off
    shading(axesArray(var), 'interp');
    axis(axesArray(var), 'equal', 'tight');
    xlabel(axesArray(var), '$x$', 'Interpreter','latex');
    ylabel(axesArray(var), '$y$', 'Interpreter','latex');
    xlim(axesArray(var),[0 L])
    ylim(axesArray(var),[0 H])
    colorbar(axesArray(var));
    titles(var) = title(axesArray(var), ...
    sprintf('M = %d', M),'Interpreter','latex');
end
drawnow;