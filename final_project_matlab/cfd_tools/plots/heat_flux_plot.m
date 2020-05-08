function fig = heat_flux_plot(heat_flux, time, filename)

fig = figure(2);
fig.Units = "inches";
fig.Position = [5, 1.5, 10, 5];
clf(fig,'reset');

subplot(1, 2, 1);
plot(time, heat_flux.left);
title("Left Wall Heat Flux");
xlabel("Time");
ylabel( ...
    "$\frac{Q}{A}=\frac{\sum_{j}^{N}\theta_{0.5,j} - \theta_{0,j}}{N}$", ...
    'Interpreter','latex', ...
    'FontSize', 15 ...
);

subplot(1, 2, 2);
plot(time, heat_flux.right);
title("Right Wall Heat Flux");
xlabel("Time");
ylabel( ...
    "$\frac{Q}{A}=\frac{\sum_{j}^{N}\theta_{M,j} - \theta_{M-0.5,j}}{N}$", ...
    'Interpreter','latex', ...
    'FontSize', 15 ...
);

saveas(fig, filename)

end
