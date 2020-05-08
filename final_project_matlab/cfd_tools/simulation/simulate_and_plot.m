function [u_vel, v_vel, theta, heat_flux, time] = simulate_and_plot( ...
    x_len, y_len, x_num, y_num, dx, dy, k_mod, ghost, ...
    u_vel, v_vel, theta, ...
    left_wall_temperature, right_wall_temperature, heat_flux, ...
    cfl_target, time, start_time, end_time, ...
    number_of_plots, current_plot_number, ...
    prandtl, rayleigh, ...
    output_dir, plot_final_state ...
)

%%%%%%%%%%%%%%%%%%%% Additional setup %%%%%%%%%%%%%%%%%%%%%%
% `plot_final_state` is an optional flag indicating to plot the final state.
% Default will generate a final plot
try
    plot_final_state;
catch
    plot_final_state = 1;
end
% Create output directory
mkdir(output_dir)

%find boundaries
gp1 = ghost + 1;
top_wall = ghost + y_num;
right_wall = ghost + x_num;

% X and Y arrays
u_x = 0:dx:x_len;
u_y = dy/2:dy:y_len-dy/2;

v_x = dx/2:dx:x_len-dx/2;
v_y = 0:dy:y_len;

theta_x = v_x;
theta_y = u_y;


% Calculate initial divergence
vel_divergence = velocity_divergence( ...
    u_vel, v_vel, dx, dy, x_num, y_num, ghost ...
);

% Calculate initial kinetic energy
kinetic_energy = get_kinetic_energy( ...
    u_vel, v_vel, top_wall, right_wall, ghost ...
);

% Calculate initial thermal energy
thermal_energy = get_thermal_energy(theta);
thermal_energy = thermal_energy(gp1:right_wall, gp1:top_wall);

% plot initial condition
file_name = output_dir + "/" + num2str(current_plot_number) ...
    + "_initial_condition.png";
simulation_plot( ...
    u_vel, u_x, u_y, ...
    v_vel, v_x, v_y, ...
    theta(gp1:right_wall, gp1:top_wall),  theta_x, theta_y,...
    0, ...time
    x_num, y_num, ghost, ...
    vel_divergence, kinetic_energy, thermal_energy, ...
    file_name ...
);

iterations = number_of_plots - 1;
time_chunk = end_time / iterations;
current_time = start_time;

%%%%%%%%%%%%%%%%%%%% Advance Simulation %%%%%%%%%%%%%%%%%%%%%%

for n=1:iterations

    chunk_end_time = n * time_chunk;

    [u_vel, v_vel, theta, heat_flux, time] = simulate( ...
        u_vel, v_vel, theta, heat_flux, dx, dy, ...
        0, cfl_target, time, current_time, chunk_end_time, ...
        prandtl, rayleigh, ...
        x_num, y_num, left_wall_temperature, right_wall_temperature,...
        ghost, k_mod  ...
    );

    % Update current time
    current_time = max(time);

    % Calculate initial divergence
    vel_divergence = velocity_divergence( ...
        u_vel, v_vel, dx, dy, x_num, y_num, ghost ...
    );

    % Calculate initial kinetic energy
    kinetic_energy = get_kinetic_energy( ...
        u_vel, v_vel, top_wall, right_wall, ghost ...
    );

    % Calculate initial thermal energy
    thermal_energy = get_thermal_energy(theta);
    thermal_energy = thermal_energy(gp1:right_wall, gp1:top_wall);

    % Output file name
    plot_num = current_plot_number + n;
    filename = num2str(plot_num) + "_time_" + num2str(chunk_end_time, "%.3g");

    simulation_plot( ...
        u_vel, u_x, u_y, ...
        v_vel, v_x, v_y, ...
        theta(gp1:right_wall, gp1:top_wall),  theta_x, theta_y,...
        current_time, ...
        x_num, y_num, ghost, ...
        vel_divergence, kinetic_energy, thermal_energy, ...
        output_dir + "/" + replace(filename, ".", "p") + ".png" ...
    );

    fprintf("Progress : %.2d \n", 100*n/iterations)

end

heat_flux_plot(heat_flux, time, output_dir + "/heat_flux.png");

if plot_final_state

    filename = num2str(number_of_plots) + "_time_" + num2str(end_time, "%.3g");
    simulation_plot( ...
        u_vel, u_x, u_y, ...
        v_vel, v_x, v_y, ...
        theta(gp1:right_wall, gp1:top_wall),  theta_x, theta_y,...
        current_time, ...
        x_num, y_num, ghost, ...
        vel_divergence, kinetic_energy, thermal_energy, ...
        output_dir + "/" + replace(filename, ".", "p") + ".png" ...
    );

end

end
