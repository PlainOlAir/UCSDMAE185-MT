%% --- Animator ---
tile = tiledlayout(3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% initialize handles
axesArray = gobjects(1, 7);
h = gobjects(1, 6);
convergenceh = gobjects(1, 4);
var_labels = {'$\rho$', '$u$', '$v$', '$e$', '$p$', '$T$', '$Convergence$'};
titles = gobjects(1, 7);
convergence_vals = output_vars{7}; 
time = time(:);

% initial plot for each variable
for var = 1:7
    if var == 7
        axesArray(var) = nexttile(tile, var, [1 3]);
        % For convergence plot (line plot)
        for k = 1:4
            convergenceh(k) = plot(axesArray(var), convergence_vals(1,k));
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
        % For field variables (pcolor plots)
        output_frame = squeeze(output_vars{var}(:,:,1));
        h(var) = pcolor(axesArray(var), xx, yy, output_frame);
        shading(axesArray(var), 'interp');
        axis(axesArray(var), 'equal', 'tight');
        xlabel(axesArray(var), '$x$', 'Interpreter','latex');
        ylabel(axesArray(var), '$y$', 'Interpreter','latex');
        colorbar(axesArray(var));
        titles(var) = title(var_labels{var}, 'Interpreter', 'latex');
    end
    
    tiletitle = title(tile, ...
        sprintf('$t=%.4e$ s\n(%d/%d)', time(1), 1, step_total),'Interpreter','latex');
end

for i = 1:10:step_total
    for var = 1:7
        if var == 7
            for k = 1:4
                % Update convergence plot using refreshdata
                refreshdata(convergenceh(k), 'caller');
                % Adjust axes limits if needed
                if i > 1
                    axesArray(var).XLim = [0, i];
                    ylims = [min(convergence_vals(1:i,:), [], 'all'), max(convergence_vals(1:i, :), [], 'all')];
                    if ylims(1) ~= ylims(2)  % Only adjust if not constant
                        axesArray(var).YLim = ylims;
                    end
                end
            end
        else
            % Update field plots directly
            output_frame = squeeze(output_vars{var}(:,:,i));
            h(var).CData = output_frame;
        end
        % Update titles
        tiletitle.String = sprintf('$t=%.3f$ ns (%d/%d)', time(i) * 1e9, i, step_total);
    end
    drawnow;
end