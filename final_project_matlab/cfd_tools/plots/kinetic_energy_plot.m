function fig = kinetic_energy_plot(kinetic_trend, initial_ke, filename)

fig = figure("Units", "inches", "Position", [0, 0, 14, 8]);

initial = [min(kinetic_trend(1, :)), max(kinetic_trend(1, :))];

subplot(2,1,1)
plot(kinetic_trend(1, :), kinetic_trend(3, :))
hold on
plot(initial, [initial_ke, initial_ke])
title("Total Kinetic Energy")
xlabel("dt")
ylabel("Kinetic Energy")
legend("Simulation KE", "Initial KE")

log_fit_type = fittype( ...
    'a + b*x.^k', ...
    'dependent', {'y'}, ...
    'independent', {'x'}, ...
    'coefficients', {'a', 'b', 'k'} ...
);

convergence_fit = fit( ...
    kinetic_trend(1, :)', ...
    kinetic_trend(2, :)', ...
    log_fit_type,  ...
    'StartPoint', [min(kinetic_trend(2, :)), 1, 3] ...
);

conv_fit = convergence_fit(kinetic_trend(1, :));

subplot(2,1,2)
loglog(kinetic_trend(1, :), kinetic_trend(2, :), 'xr')
hold on
loglog(kinetic_trend(1, :), conv_fit)
str_title = sprintf( ...
    "Kinetic Energy Error: Convergence = %.2g", convergence_fit.k ...
);
title(str_title)
xlabel("dt")
ylabel("Kinetic Energy Loss")
legend("KE Loss", "Convergence Fit")

saveas(fig, filename)

end

