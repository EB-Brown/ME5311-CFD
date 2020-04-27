function fig = energy_plot( ...
    kinetic_trend, initial_ke, ...
    thermal_trend, initial_thermal, ...
    filename ...
)

fig = figure("Units", "inches", "Position", [0, 0, 14, 8]);

dt_range = [min(kinetic_trend(1, :)), max(kinetic_trend(1, :))];

%%%%%%%%%%%%%%%%%%%%%%%%%%% Total Energy Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% Kinetic Energy %%%%%%%
ax = subplot(2,2,1);
ax.Position = [0.13,0.525,.33,.425];

plot(kinetic_trend(1, :), kinetic_trend(3, :))
hold on
plot(dt_range, [initial_ke, initial_ke])

min_ke = min(kinetic_trend(3, :));
ke_plot_tol = (initial_ke - min_ke) * 0.1;
ylim([min_ke - ke_plot_tol, initial_ke + ke_plot_tol])

title("Total Kinetic Energy")
xlabel("dt")
ylabel("Kinetic Energy")
legend("Simulation KE", "Initial KE", 'Location', 'southoutside')

%%%%%%% Thermal Energy %%%%%%%
ax = subplot(2,2,2);
ax.Position = [0.58,0.525,.33,.425];

plot(thermal_trend(1, :), thermal_trend(3, :))
hold on
plot(dt_range, [initial_thermal, initial_thermal])

min_therm = min(thermal_trend(3, :));
therm_plot_tol = (initial_thermal - min_therm) * 0.1;
ylim([min_therm - therm_plot_tol, initial_thermal + therm_plot_tol])

title("Total Thermal Energy")
xlabel("dt")
ylabel("Thermal Energy")
legend( ...
    "Simulation Thermal Energy", ...
    "Initial Thermal Energy", ...
    'Location', ...
    'southoutside' ...
)

%%%%%%%%%%%%%%%%%%%%%%%%%%% Convergence fitting %%%%%%%%%%%%%%%%%%%%%%%%%%%

log_fit_type = fittype( ...
    'a + b*x.^k', ...
    'dependent', {'y'}, ...
    'independent', {'x'}, ...
    'coefficients', {'a', 'b', 'k'} ...
);

kinetic_fit = fit( ...
    kinetic_trend(1, :)', ...
    kinetic_trend(2, :)', ...
    log_fit_type,  ...
    'StartPoint', [min(kinetic_trend(2, :)), 6.659e+07, 3] ...
);

thermal_fit = fit( ...
    thermal_trend(1, :)', ...
    thermal_trend(2, :)', ...
    log_fit_type,  ...
    'StartPoint', [min(thermal_trend(2, :)), 3.884e+08, 3] ...
);

%%%%%%%%%%%%%%%%%%%%%%%%%%% Convergence Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,2,3)
loglog(kinetic_trend(1, :), kinetic_trend(2, :), 'xr')

kinetic_title = sprintf( ...
    "Kinetic Energy Error:\nConvergence = %.2g", kinetic_fit.k ...
);
title(kinetic_title)
xlabel("dt")
ylabel("Kinetic Energy Loss")

subplot(2,2,4)
loglog(thermal_trend(1, :), thermal_trend(2, :), 'xr')

thermal_title = sprintf( ...
    "Thermal Energy Error\nConvergence = %.2g", thermal_fit.k ...
);
title(thermal_title)
xlabel("dt")
ylabel("Thermal Energy Loss")

%%%%%%%%%%%%%%%%%%%%%%%%%%% Save Figure %%%%%%%%%%%%%%%%%%%%%%%%%%%
saveas(fig, filename)

end

