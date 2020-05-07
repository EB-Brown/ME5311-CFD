clear all
clc

%%%%%%%%%%%%%%%%% Inputs %%%%%%%%%%%%%%%%%

% physical domain
x_len = 0.5;
y_len = 2;
end_time = 0.015;
left_wall_temperature = 1; % const or func where 0 <= left_wall_temperature <= 1

% Dynamic time steps are calculated using the relaxation method
cfl_target = 0.7;

% Number of grid cells
x_num = 64;
y_num = 256;
ghost = 1; % number  of ghost cells

% Simulation Parameters
prandtl = 0.7; % Prandtl number for air
rayleigh = 1e6; % Rayleigh number controls convective driver

% Number of plots generated in equal time intervals.
% Initial and final profiles are included
number_of_plots = 11;

%%%%%%%%%%%%%%%%% Setup %%%%%%%%%%%%%%%%%

% find boundaries
gp1 = ghost + 1;
top_wall = ghost + y_num;
right_wall = ghost + x_num;

% Calculate space steps
dx = x_len/x_num;
dy = y_len/y_num;

% X and Y arrays
u_x = 0:dx:x_len;
u_y = dy/2:dy:y_len-dy/2;

v_x = dx/2:dx:x_len-dx/2;
v_y = 0:dy:y_len;

theta_x = v_x;
theta_y = u_y;

% Calculate matrix for poisson solution
k_mod = get_k_mod(x_num, y_num, dx, dy);

% Fill u, v, and temperature with random numbers
u_velocity = zeros(x_num + 2 * ghost, y_num + 2 * ghost);
v_velocity = zeros(x_num + 2 * ghost, y_num + 2 * ghost);

temperature = zeros(x_num + 2 * ghost, y_num + 2 * ghost) + 0.5;

% Set velocity boundary condition
[u_velocity, v_velocity] = velocity_boundaries( ...
    u_velocity, v_velocity, x_num, y_num, ghost ...
);

% Set temperature boundary condition
temperature = temperature_boundaries( ...
    temperature, x_num, y_num, ghost, left_wall_temperature ...
);

% Remove Divergence
[u_velocity, v_velocity] = remove_divergence( ...
    u_velocity, v_velocity, dx, dy, x_num, y_num, k_mod, ghost ...
);

% Calculate initial divergence
vel_divergence = velocity_divergence( ...
    u_velocity, v_velocity, dx, dy, x_num, y_num, ghost ...
);

fprintf( ...
    "max initial divergence = %.2g \n", ...
    max(max(abs(vel_divergence))) ...
);

% Calculate initial kinetic energy
kinetic_energy = get_kinetic_energy( ...
    u_velocity, v_velocity, top_wall, right_wall, ghost ...
);

% Calculate initial thermal energy
thermal_energy = get_thermal_energy(temperature);
thermal_energy = thermal_energy(gp1:right_wall, gp1:top_wall);
%%%%%%%%%%%%%%%%% Simulation %%%%%%%%%%%%%%%%%

output_dir = "plots/simple_boundary_conditions";
scaled_parms = "/"+ "pr_" + num2str(prandtl, "%.2e") ...
    + "__ra_" + num2str(rayleigh, "%.2e");
output_dir = output_dir + replace(replace(scaled_parms, ".", "_"), "+", "");
mkdir(output_dir)

simulation_plot( ...
    u_velocity, u_x, u_y, ...
    v_velocity, v_x, v_y, ...
    temperature(gp1:right_wall, gp1:top_wall),  theta_x, theta_y,...
    0, ...
    x_num, y_num, ghost, ...
    vel_divergence, kinetic_energy, thermal_energy, ...
    output_dir + "/" + "0_initial_condition.png" ...
);

iterations = number_of_plots - 1;
time_chunk = end_time / iterations;

time = 0;
for n=1:iterations

    chunk_end_time = n * time_chunk;

    [u_velocity, v_velocity, temperature, time] = simulate( ...
        u_velocity, v_velocity, temperature, dx, dy, ...
        0, cfl_target, time, chunk_end_time, ...
        prandtl, rayleigh, ...
        x_num, y_num, left_wall_temperature, ...
        ghost, k_mod  ...
    );

    % Calculate initial divergence
    vel_divergence = velocity_divergence( ...
        u_velocity, v_velocity, dx, dy, x_num, y_num, ghost ...
    );

    fprintf( ...
        "max divergence = %.2g \n", ...
        max(max(abs(vel_divergence))) ...
    );

    % Calculate initial kinetic energy
    kinetic_energy = get_kinetic_energy( ...
        u_velocity, v_velocity, top_wall, right_wall, ghost ...
    );

    % Calculate initial thermal energy
    thermal_energy = get_thermal_energy(temperature);
    thermal_energy = thermal_energy(gp1:right_wall, gp1:top_wall);

    filename = num2str(n) + "_time_" + num2str(chunk_end_time, "%.3g");

    simulation_plot( ...
        u_velocity, u_x, u_y, ...
        v_velocity, v_x, v_y, ...
        temperature(gp1:right_wall, gp1:top_wall),  theta_x, theta_y,...
        time, ...
        x_num, y_num, ghost, ...
        vel_divergence, kinetic_energy, thermal_energy, ...
        output_dir + "/" + replace(filename, ".", "p") + ".png" ...
    );

    fprintf("Progress : %.2d \n", 100*n/iterations)

end

