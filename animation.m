%% --- Animator ---
tile = tiledlayout(3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% initialize handles
axesArray = gobjects(1, 7);
h = gobjects(1, 6);
convergenceh = gobjects(1, 4);
var_labels = {'$\rho$', '$u$', '$v$', '$e$', '$p$', '$T$', '$Convergence$'};
titles = gobjects(1, 7);
convergence_vals = output_vars{7}; 

% initial plot for each variable
for var = 1:7
    if var == 7
        axesArray(var) = nexttile(tile, var, [1 3]);
        for k = 1:4
            convergenceh(k) = semilogy(axesArray(var), convergence_vals(1,k));
            set(convergenceh(k), 'XDataSource', '1:i')
            set(convergenceh(k), 'YDataSource', 'convergence_vals(1:i,k)');
            hold on
        end
        hold off
        xlabel(axesArray(var), '$step$', 'Interpreter','latex');
        ylabel(axesArray(var), '$Convergence$', 'Interpreter','latex');
        legend({'$\rho$', '$\rho u$', '$\rho v$', '$E_t$'}, 'Interpreter', 'latex')
    else
        axesArray(var) = nexttile(tile, var);
        output_frame = squeeze(output_vars{var}(:,:,1));
        h(var) = pcolor(axesArray(var), xx, yy, output_frame);
        shading(axesArray(var), 'interp');
        axis(axesArray(var), 'equal', 'tight');
        xlabel(axesArray(var), '$x$', 'Interpreter','latex');
        ylabel(axesArray(var), '$y$', 'Interpreter','latex');
        colorbar(axesArray(var));
        titles(var) = title(var_labels{var}, 'Interpreter', 'latex');
    end
    tiletitle = title(tile, sprintf('Step %d/%d', i, step_total));
end
% draw plots at time stepping 50
for i = 1:100:step_total
    for var = 1:7
        if var == 7
            for k = 1:4
                refreshdata(convergenceh(k), 'caller');
                if i > 1
                    axesArray(var).XLim = [0, i];
                    ylims = [min(convergence_vals(1:i,:), [], 'all'), max(convergence_vals(1:i, :), [], 'all')];
                    if ylims(1) ~= ylims(2)  
                        axesArray(var).YLim = ylims;
                    end
                end
            end
        else
            output_frame = squeeze(output_vars{var}(:,:,i));
            h(var).CData = output_frame;
        end
        tiletitle.String = sprintf('Step %d/%d', i, step_total);
    end
    drawnow;
end
for var = 1:7
    if var == 7
        axesArray(var) = nexttile(tile, var, [1 3]);
for k = 1:4
    convergenceh(k) = semilogy(axesArray(var), 1:1500, convergence_vals(1:1500,k));
    set(convergenceh(k), 'XDataSource', '1:1500');
    set(convergenceh(k), 'YDataSource', sprintf('convergence_vals(1:1500,%d)', k));
    hold on
end
        hold off
        xlabel(axesArray(var), '$step$', 'Interpreter','latex');
        ylabel(axesArray(var), '$Convergence$', 'Interpreter','latex');
        legend({'$\rho$', '$\rho u$', '$\rho v$', '$E_t$'}, 'Interpreter', 'latex')
    else
        axesArray(var) = nexttile(tile, var);
        output_frame = squeeze(output_vars{var}(:,:,end-1));
        h(var) = pcolor(axesArray(var), xx, yy, output_frame);
        shading(axesArray(var), 'interp');
        axis(axesArray(var), 'equal', 'tight');
        xlabel(axesArray(var), '$x$', 'Interpreter','latex');
        ylabel(axesArray(var), '$y$', 'Interpreter','latex');
        colorbar(axesArray(var));
        titles(var) = title(var_labels{var}, 'Interpreter', 'latex');
    end
    tiletitle = title(tile, sprintf('Step %d/%d', step_total, step_total));
end
drawnow;