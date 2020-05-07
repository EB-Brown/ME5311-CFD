function fig = simulation_plot( ...
    u_velocity, u_x, u_y, ...
    v_velocity, v_x, v_y, ...
    temperature,  theta_x, theta_y,...
    time, ...
    x_num, y_num, ghost, ...
    vel_divergence, kinetic_energy, thermal_energy, filename ...
)

gp1 = ghost + 1;
roof = ghost + y_num;
right = ghost + x_num;

x = 0.21875;
y = 0.35;
x_size = x;
y_size = y;
x_shift = 0;
y_shift = 0;

aspect_ratio = x_size / y_size;

domain_ratio = max(u_x) / max(u_y);
if domain_ratio > aspect_ratio
    y_size = x_size / domain_ratio;
    y_shift = (0.45-0.05) * ((y-y_size)/(2*y));
elseif domain_ratio < aspect_ratio
    x_size = y_size * domain_ratio;
    x_shift = (0.33-0.065625) * ((x-x_size)/(2*x));
end

fig = figure(1);
fig.Units = "inches";
fig.Position = [0, 0, 14, 8];
clf(fig,'reset');
annotation( ...
    'textbox', [0.3954,0.975,0.21,0], ...
    'string', sprintf("time = %g", time), ...
    'LineStyle', 'none', ...
    'FontSize', 16, ...
    'FontWeight', 'bold' ...
);

ax = subplot(2,3,1);
ax.Position = [0.065625+x_shift,0.525+y_shift,x_size,y_size];
contourf( ...
    u_x, u_y, u_velocity(ghost:right, gp1: roof)', ...
    100, ...
    'LineStyle', 'none' ...
)
colormap(autumn)
title("U Velocity")
colorbar()


ax = subplot(2,3,4);
ax.Position = [0.065625+x_shift,0.05+y_shift,x_size,y_size];
contourf( ...
    v_x, v_y, v_velocity(gp1:right, ghost: roof)', ...
    100, ...
    'LineStyle', 'none' ...
)
colormap(autumn)
title("V Velocity")
colorbar()

ax = subplot(2,3,2);
ax.Position = [0.38+x_shift,0.525+y_shift,x_size,y_size];
contourf(theta_x, theta_y, vel_divergence', 100, 'LineStyle', 'none')
colormap(autumn)
divergence_title = sprintf( ...
    "Velocity Divergence\nMax Abs = %.2e",max(max(abs(vel_divergence))) ...
);
title(divergence_title)
colorbar()

ax = subplot(2,3,3);
ax.Position = [0.71+x_shift,0.525+y_shift,x_size,y_size];
contourf(u_x, v_y, kinetic_energy', 100, 'LineStyle', 'none')
colormap(autumn)
title("Kinetic Energy")
colorbar()

ax = subplot(2,3,5);
ax.Position = [0.38+x_shift,0.05+y_shift,x_size,y_size];
contourf(theta_x, theta_y, temperature', 100, 'LineStyle', 'none')
caxis([0,1])
colormap(autumn)
title("Theta Temperature")
colorbar()

ax = subplot(2,3,6);
ax.Position = [0.71+x_shift,0.05+y_shift,x_size,y_size];
contourf(theta_x, theta_y, thermal_energy', 100, 'LineStyle', 'none')
colormap(autumn)
title("Thermal Energy")
colorbar()

saveas(fig, filename)

end

