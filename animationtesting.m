%% --- Animator ---
tile = tiledlayout(2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
% initialize handles
axesArray = gobjects(1, 8);
h = gobjects(1, 8);
var_labels = {'$\rho$', '$u$', '$v$', '$e$', '$p$', '$T$', '$Convergence$'};
titles = gobjects(1, 8);
convergence_vals = abs(squeeze(output_vars{7}(:))); 
time = time(:);

% initial plot for each variable
for var = 1:7
    axesArray(var) = nexttile(tile, var);
    if var == 7
        % For convergence plot (line plot)
        h(var) = plot(axesArray(var), time(1), convergence_vals(1), 'b-');
        xlabel(axesArray(var), '$t$', 'Interpreter','latex');
        ylabel(axesArray(var), '$Convergence$', 'Interpreter','latex');
        set(h(var), 'XDataSource', 'time(1:i)');
        set(h(var), 'YDataSource', 'convergence_vals(1:i)');
    else
        % For field variables (pcolor plots)
        output_frame = squeeze(output_vars{var}(:,:,1));
        h(var) = pcolor(axesArray(var), xx, yy, output_frame);
        shading(axesArray(var), 'interp');
        axis(axesArray(var), 'equal', 'tight');
        xlabel(axesArray(var), '$x$', 'Interpreter','latex');
        ylabel(axesArray(var), '$y$', 'Interpreter','latex');
        colorbar(axesArray(var));
    end
    
    titles(var) = title(axesArray(var), ...
        sprintf('%s at $t=%.4e$ s\n(%d/%d)', var_labels{var}, time(1), 1, step_total),'Interpreter','latex');
end

for i = 1:50:step_total
    for var = 1:7
        if var == 7
            % Update convergence plot using refreshdata
            refreshdata(h(var), 'caller');
            % Adjust axes limits if needed
            if i > 1
                axesArray(var).XLim = [time(1), time(i)];
                ylims = [min(convergence_vals(1:i)), max(convergence_vals(1:i))];
                if ylims(1) ~= ylims(2)  % Only adjust if not constant
                    axesArray(var).YLim = ylims;
                end
            end
        else
            % Update field plots directly
            output_frame = squeeze(output_vars{var}(:,:,i));
            h(var).CData = output_frame;
        end
        % Update titles
        titles(var).String = sprintf('%s at $t=%.4e$ s\n(%d/%d)', var_labels{var}, time(i), i, step_total);
    end
    drawnow;
end