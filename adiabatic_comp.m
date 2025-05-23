% Implement the adiabatic wall boundary condition into your code. Make sure that 2nd-order accuracy is maintained. 
% Run the Mach 4 case for an adiabatic wall and compare the normalized temperature and pressure profiles (as a function of y)
% at three different x-locations, x/L = 0.25, 0.5, 0.75.
% Compare the wall temperature (as a function of x) for both cases in a second plot.

xloc = [0.25, 0.5, 0.75];
x_indices = round(xloc * nx);

output_vars_adiabatic = struct2cell(load("adiabatic_output.mat"));
output_vars_adiabatic = [output_vars_adiabatic{:}];
output_vars = struct2cell(load("output.mat"));
output_vars = [output_vars{:}];
tile = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile(1)
hold on
for x_iter = 1:3
    plot(linspace(0,H,ny), output_vars{5}(x_indices(1,x_iter),:,1500)/pinf);
    plot(linspace(0,H,ny), output_vars_adiabatic{5}(x_indices(1,x_iter),:,1500)/pinf,'LineStyle','--');
end
xlabel('$y$', 'Interpreter','latex');
ylabel('$p_0$', 'Interpreter','latex');
xlim([0 H])
title(sprintf('Pressure at $x$ = 0.25, 0.50, 0.75'),'Interpreter','latex');
legend('$T_{inf}$ $x$ = 0.25','Adiabatic $x$ = 0.25','$T_{inf}$ $x$ = 0.50','Adiabatic $x$ = 0.50','$T_{inf}$ = 0.75','Adiabatic $x$ = 0.75','Interpreter','latex')

nexttile(2)
hold on
for x_iter = 1:3
    plot(linspace(0,H,ny), output_vars{6}(x_indices(1,x_iter),:,1500)/Tinf);
    plot(linspace(0,H,ny), output_vars_adiabatic{6}(x_indices(1,x_iter),:,1500)/Tinf,'LineStyle','--');
end
xlabel('$y$', 'Interpreter','latex');
ylabel('$T_0$', 'Interpreter','latex');
xlim([0 H])
title(sprintf('Temperature at $x$ = 0.25, 0.50, 0.75'),'Interpreter','latex');
legend('$p_{inf}$ $x$ = 0.25','Adiabatic $x$ = 0.25','$p_{inf}$ $x$ = 0.50','Adiabatic $x$ = 0.50','$p_{inf}$ $x$ = 0.75','Adiabatic $x$ = 0.75','Interpreter','latex')

nexttile(3)
hold on
plot(linspace(0,L,nx), output_vars{6}(:,1,1500));
plot(linspace(0,L,nx), output_vars_adiabatic{6}(:,1,1500),'LineStyle','--');
xlabel('$x$', 'Interpreter','latex');
ylabel('$Wall T$', 'Interpreter','latex');
xlim([0 H])
title(sprintf('Wall Temperature'),'Interpreter','latex');
legend('$T_{inf}$ $T$','Adiabatic $T$','Interpreter','latex')