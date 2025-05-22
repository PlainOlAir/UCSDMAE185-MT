%% --- Animator ---
tile = tiledlayout(2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
% initialize handles
axesArray = gobjects(1, 8);
h = gobjects(1, 8);
var_labels = {'$\rho$', '$u$', '$v$', '$e$', '$p$', '$T$', '$Convergence$'};
titles = gobjects(1, 8);

% initial plot for each variable
for var = 1:12
    if var >= 7
        % For convergence plots (line plots) - create just one axes for all convergence lines
        if var == 7
            axesArray(7) = nexttile(tile, 7);
            hold(axesArray(7), 'on');
            xlabel(axesArray(7), '$t$', 'Interpreter','latex');
            ylabel(axesArray(7), '$log_{10}(Residual)$', 'Interpreter','latex');
            titles(7) = title(axesArray(7), ...
                sprintf('Convergence Metrics at $t=%.4e$ s\n(%d/%d)', time(1), 1, step_total),'Interpreter','latex');
            
            % Initialize all 6 convergence lines
            colors = lines(6); % Get 6 distinct colors
            h_conv = gobjects(6,1); % Preallocate graphics objects
            
            for conv_var = 1:6
                output_data = squeeze(output_vars{6+conv_var}(:,:,1));
                h_conv(conv_var) = semilogy(axesArray(7), time(1), output_data, '-', ...
                    'Color', colors(conv_var,:), 'DisplayName', var_labels{conv_var});
            end
            legend(axesArray(7), 'show', 'Interpreter', 'latex');
            set(axesArray(7), 'YScale', 'log'); % Ensure y-axis is logarithmic
            grid(axesArray(7), 'on');
        end
    else
        % For field variables (pcolor plots)
        axesArray(var) = nexttile(tile, var);
        output_frame = squeeze(output_vars{var}(:,:,1));
        h(var) = pcolor(axesArray(var), xx, yy, output_frame);
        shading(axesArray(var), 'interp');
        axis(axesArray(var), 'equal', 'tight');
        xlabel(axesArray(var), '$x$', 'Interpreter','latex');
        ylabel(axesArray(var), '$y$', 'Interpreter','latex');
        colorbar(axesArray(var));
        titles(var) = title(axesArray(var), ...
            sprintf('%s at $t=%.4e$ s\n(%d/%d)', var_labels{var}, time(1), 1, step_total),'Interpreter','latex');
    end
end

for i = 1:50:step_total
    for var = 1:12
        if var >= 7
            % Update convergence plot title only once
            if var == 7
                titles(7).String = sprintf('Convergence Metrics at $t=%.4e$ s\n(%d/%d)', time(i), i, step_total);
                % Update all convergence lines
                for conv_var = 1:6
                    output_data = squeeze(output_vars{6+conv_var}(:,:,1:i));
                    set(h_conv(conv_var), 'XData', time(1:i), 'YData', output_data);
                end
            end
        else
            % Update field plots directly
            output_frame = squeeze(output_vars{var}(:,:,i));
            h(var).CData = output_frame;
            titles(var).String = sprintf('%s at $t=%.4e$ s\n(%d/%d)', var_labels{var}, time(i), i, step_total);
        end
    end
    drawnow;
end