%% --- Animator ---
tile = tiledlayout(3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% initial plot for each variable
for var = 1:6
    nexttile(var)
    pcolor(xx, yy, squeeze(output_vars{var}(:, :, 1)))
    shading interp
    axis equal tight
end

nexttile(7, [1 3])
convPlot = plot(output_vars{7});

for i = 1:10:step_total
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