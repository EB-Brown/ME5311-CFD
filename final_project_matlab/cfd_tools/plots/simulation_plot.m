function fig = simulation_plot( ...
    u_velocity, u_x, u_y, ...
    v_velocity, v_x, v_y, ...
    temperature,  theta_x, theta_y,...
    x_num, y_num, ghost, ...
    vel_divergence, kinetic_energy, thermal_energy, filename ...
)

fig = figure("Units", "inches", "Position", [0, 0, 14, 8]);

subplot(2,3,1)
contour( ...
    u_x, u_y, u_velocity(ghost:ghost + x_num, ghost + 1: ghost + y_num)', 100 ...
)
title("U Velocity")
colorbar()

subplot(2,3,4)
contour( ...
    v_x, v_y, v_velocity(ghost + 1:ghost + x_num, ghost: ghost + y_num)', 100 ...
)
title("V Velocity")
colorbar()

subplot(2,3,2)
contour(theta_x, theta_y, vel_divergence', 100)
divergence_title = sprintf( ...
    "Velocity Divergence\nMax Abs = %.2e",max(max(abs(vel_divergence))) ...
);
title(divergence_title)
colorbar()

subplot(2,3,3)
contour(u_x, v_y, kinetic_energy', 100)
title("Kinetic Energy")
colorbar()

subplot(2,3,5)
contour(theta_x, theta_y, temperature', 100)
title("Theta Temperature")
colorbar()

subplot(2,3,6)
contour(theta_x, theta_y, thermal_energy', 100)
title("Thermal Energy")
colorbar()

saveas(fig, filename)

end

