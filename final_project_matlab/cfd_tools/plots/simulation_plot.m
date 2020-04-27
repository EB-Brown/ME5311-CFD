function fig = simulation_plot( ...
    u_velocity, v_velocity, x_num, y_num, ghost, ...
    vel_divergence, kinetic_energy, filename ...
)

fig = figure("Units", "inches", "Position", [0, 0, 14, 8]);

subplot(2,2,1)
contour(u_velocity(ghost:ghost + x_num, ghost + 1: ghost + y_num), 100)
title("U Velocity")
colorbar()

subplot(2,2,3)
contour(v_velocity(ghost + 1:ghost + x_num, ghost: ghost + y_num), 100)
title("V Velocity")
colorbar()

subplot(2,2,2)
contour(vel_divergence, 100)
title("Velocity Divergence")
colorbar()

subplot(2,2,4)
contour(kinetic_energy, 100)
title("Kinetic Energy")
colorbar()

saveas(fig, filename)

end

