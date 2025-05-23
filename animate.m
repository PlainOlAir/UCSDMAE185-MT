%%%%%%%%%%%%%%%%
%%% Animator %%%
%%%%%%%%%%%%%%%%

%% Setup

% Open figure
tile = tiledlayout(3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% Resize figure
set(gcf, "Position", [100, 100, animationwidth, animationheight]);

% Cell array for plot labels
plot_labels = {'$\rho ~ (kg/m^3)$', '$u ~ (m/s)$', '$v ~ (m/s)$', '$e ~ (J/kg)$', '$p ~ (Pa)$', '$T ~ (^\circ K)$', 'Normalized Residuals'};

%% Preallocation

% Allocate handle arrays for graphical objects
axesArray = gobjects(1, 7);
h = gobjects(1, 6);
residualh = gobjects(1, 4);
titles = gobjects(1, 7);
animationsteps = [1 animationstep+1:animationstep:step_total+1 step_total+1];
im = cell(1,length(animationsteps));

%% Initialization

% Initialize plots
for var = 1:7

    % Residual plot
    if var == 7

        % Span residual plot across the bottom
        axesArray(var) = nexttile(tile, var, [1 3]);

        % Plot residuals of all 4 conservative variables on a logarithmic scale
        for k = 1:4
            residualh(k) = semilogy(axesArray(var), output_vars{7}(1,k));

            % Setting data source to easily call refreshdata
            set(residualh(k), 'XDataSource', '1:animationsteps(i)');
            set(residualh(k), 'YDataSource', 'output_vars{7}(1:animationsteps(i),k)');
            hold on
        end
        hold off

        % Axes label & legend
        xlabel(axesArray(var), 'step');
        ylabel(axesArray(var), 'Residual');
        legend({'$\rho$', '$\rho u$', '$\rho v$', '$E_t$'}, 'Interpreter', 'latex')
        title(axesArray(var), plot_labels{var})

        % Primitive variable plots
    else
        axesArray(var) = nexttile(tile, var);
        h(var) = pcolor(axesArray(var), xx, yy, reshape(output_vars{var}(:,:,1), [nx ny]));

        % Set colormaps for plots
        colormap(axesArray(var), "turbo")
        if var == 6
            colormap(axesArray(var), "hot")
        end

        % Interpolated shading, axes label, colorbar, and titles
        shading(axesArray(var), 'interp');
        axis(axesArray(var), 'equal', 'tight');
        xlabel(axesArray(var), '$x ~ (m)$', 'Interpreter','latex');
        ylabel(axesArray(var), '$y ~ (m)$', 'Interpreter','latex');
        colorbar(axesArray(var));
        titles(var) = title(plot_labels{var}, 'Interpreter', 'latex');
    end

    % Main title for figure
    tiletitle = title(tile, sprintf('Time: %.3f ns\nStep: %d/%d',time(1) * 1e9, 1, step_total), 'Interpreter', 'latex');
end

%% Animation %%

if animationexport
    im{1} = frame2im(getframe(gcf));
end

% Step through data
for i = 2:length(animationsteps)

    % Step through all plots
    for var = 1:7

        % Residual plot
        if var == 7

            % Step through all 4 residuals
            for k = 1:4

                % Update residual plot using refreshdata
                refreshdata(residualh(k), 'caller');

                % Calculate Y axis limits
                axesArray(var).XLim = [0, animationsteps(i)];
                ylimold = [0 0];
                ylims = [min(output_vars{7}(1:animationsteps(i),:), [], 'all'), max(output_vars{7}(1:animationsteps(i), :), [], 'all')];

                % Adjust Y axis if changed
                if isequal(ylims, ylimold)
                    axesArray(var).YLim = ylims;
                    ylimold = ylims;
                end
            end
        
        % Primitive variable plots
        else
            output_frame = reshape(output_vars{var}(:,:,animationsteps(i)), [nx ny]);
            h(var).CData = output_frame;
        end

        % Update titles
        tiletitle.String = sprintf('Time: %.3f ns\nStep: %d/%d',time(animationsteps(i)) * 1e9, animationsteps(i)-1, step_total);
    end

    % Update figure
    drawnow;

    if animationexport
        im{i} = frame2im(getframe(gcf));
    end
end

% If animationexport is TRUE in 'setup.m', export to file specified there
if animationexport
    animateexport
end