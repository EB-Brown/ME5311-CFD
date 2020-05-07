clear all
clc

%{
Procedure:
This simulation initializes the temperature profile with a quarter sin wave. The
left wall has a fixed temperature of 1 and the right wall has a fixed temperatre
of 0. The sin wave is defined as sin(pi* x_array / (2 * x_length)) so that the
fluid temperature toward the left wall approaches 0 and the fluid temperature
toward the right wall approaches 1.

Hypothesis:
With the colder portion of the fluid against the hot wall and the hotter portion
of the fluid against the cold wall, natural convection will begin to develop in
the counter-clockwise direction. As the fluid along the wall changes its
temperature, the flow will begin to move in the clockwise direction. The final
profile will look similar to the `simulation_simple.m` results.

Observations:
The flow begins to move in the counter-clockwise direction. As the fluid along
the wall changes its temperature, the flow along the wall begins to move against
the initial convection direction. Vortex sheets begin to develop. The voertices
develop in multiple regions and the flow becomes chaotic. The flow settles with
most of the heat along the top of the domain with minor convection along the
walls.
%}

%%%%%%%%%%%%%%%%% Inputs %%%%%%%%%%%%%%%%%

%{
Modify the inputs of this section to control the simulation properties

PARAMETERS
x_len: the length of the domain in the x direction
y_len: the length of the domain in the y direction
end_time: the amount of time the simulation will run for
left_wall_temperature: The left wall temperature condition. An array can be
    used as long as it represents a continous function and its length is
    `y_len + 2*ghost`.
cfl_target: The target CFL for controlling the time step size
x_num: Number of grid cells in the x direction
y_num: Number of grid cells in the y direction
ghost: Number of ghost cells
prandtl: The Prandtl number for the simulation
rayleigh: The Rayleigh number for the simulation (Requires adjusting
    `cfl_target` accordingly to assure stability)
number_of_plots: The number of plots you would like to generate and save. A
    minimum of 1 is required to cover the initial condition. The time domain is
    divided into `end_time / (number_of_plots - 1)` chunks. A plot is generated
    after each chunk has been simulated.
%}

% physical domain
x_len = 1;
y_len = 2;
end_time = 0.01;
left_wall_temperature = 1;
right_wall_temperature = 0;

% Dynamic time steps are calculated using the relaxation method
cfl_target = 0.7;

% Number of grid cells
x_num = 150;
y_num = 300;
ghost = 1;

% Simulation Parameters
prandtl = 0.6;
rayleigh = 3.5*1e7;

% Number of plots generated in equal time intervals.
% Initial and final profiles are included
number_of_plots = 101;

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

temperature = zeros(x_num + 2 * ghost, y_num + 2 * ghost);
init_temp = sin(pi*theta_x / (2*x_len));
for row=gp1:top_wall
    temperature(gp1:right_wall, row) = init_temp;
end

% Set velocity boundary condition
[u_velocity, v_velocity] = velocity_boundaries( ...
    u_velocity, v_velocity, x_num, y_num, ghost ...
);

% Set temperature boundary condition
temperature = temperature_boundaries( ...
    temperature, ...
    x_num, ...
    y_num, ...
    ghost, ...
    left_wall_temperature, ...
    right_wall_temperature ...
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

output_dir = "plots/sin_temperature_conditions";
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
        x_num, y_num, left_wall_temperature, right_wall_temperature,...
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

