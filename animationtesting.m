%% --- Animator --- TODO: verify functioning
tile = tiledlayout(2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
% initialize handles
axesArray = gobjects(1, 8);
h = gobjects(1, 8);
var_labels = {'$\rho$', '$u$', '$v$', '$e$', '$p$', '$T$', '$Convergence$'};
titles = gobjects(1, 8);
% initial plot for each variable
for var = 1:7
    axesArray(var) = nexttile(tile, var);
    h(var) = pcolor(axesArray(var), xx, yy, squeeze(output_vars(var,:,:,1)));
    shading(axesArray(var), 'interp');
    axis(axesArray(var), 'equal', 'tight');
    xlabel(axesArray(var), '$x$', 'Interpreter','latex');
    ylabel(axesArray(var), '$y$', 'Interpreter','latex');
    colorbar(axesArray(var));
    titles(var) = title(axesArray(var), ...
        sprintf('%s at $t=%.11f$ s \\(%d/%d\\)', var_labels{var}, time(1), i, step_total),'Interpreter','latex');
end
% Animate over time
for i = 350
    for var = 1:7
        h(var).CData = squeeze(output_vars(var, :, :, i));
        titles(var).String = sprintf('%s at $t=%.11f$ s \\(%d/%d\\)', var_labels{var}, time(i), i, step_total);
    end
    drawnow;
end