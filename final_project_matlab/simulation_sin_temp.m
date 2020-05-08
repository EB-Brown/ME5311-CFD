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
x_len = 0.5;
y_len = 2;
end_time = 0.01;
left_wall_temperature = 1;
right_wall_temperature = 0;

% Dynamic time steps are calculated using the relaxation method
cfl_target = 0.5;

% Number of grid cells
x_num = 50;
y_num = 200;
ghost = 1;

% Simulation Parameters
prandtl = 0.7;
rayleigh = 1e6;

% Number of plots generated in equal time intervals.
% Initial and final profiles are included
number_of_plots = 11;

%%%%%%%%%%%%%%%%% Setup %%%%%%%%%%%%%%%%%

% Calculate space steps
dx = x_len/x_num;
dy = y_len/y_num;

% Calculate matrix for poisson solution
k_mod = get_k_mod(x_num, y_num, dx, dy);

% Initialize velocity as u(x,y) = v(x,y) = 0
[u_velocity, v_velocity] = initialize_velocity( ...
    x_num, y_num, dx, dy, k_mod, ghost ...
);

% Initial horizontal temperature profile
top_wall = ghost + y_num;
theta_x = dx/2:dx:x_len-dx/2;
init_temp = sin(pi*theta_x / (2*x_len));

% Assign profile to each row of temperature
gp1 = ghost + 1;
right_wall = ghost + x_num;
temperature = zeros(x_num + 2 * ghost, y_num + 2 * ghost);
for row=ghost+1:top_wall
    temperature(gp1:right_wall, row) = init_temp;
end

% Set temperature boundary condition
temperature = temperature_boundaries( ...
    temperature, ...
    x_num, ...
    y_num, ...
    ghost, ...
    left_wall_temperature, ...
    right_wall_temperature ...
);

heat_flux.left = [];
heat_flux.right = [];
time = [];

%%%%%%%%%%%%%%%%% Simulation %%%%%%%%%%%%%%%%%

output_dir = "plots/sin_temperature_conditions";
scaled_parms = "/"+ "pr_" + num2str(prandtl, "%.2e") ...
    + "__ra_" + num2str(rayleigh, "%.2e");
output_dir = output_dir + replace(replace(scaled_parms, ".", "_"), "+", "");

[u_velocity, v_velocity, theta, heat_flux, time] = simulate_and_plot( ...
    x_len, y_len, x_num, y_num, dx, dy, k_mod, ghost, ...
    u_velocity, v_velocity, temperature, ...
    left_wall_temperature, right_wall_temperature, heat_flux, ...
    cfl_target, time, 0, end_time, ...
    number_of_plots, 0, ...
    prandtl, rayleigh, ...
    output_dir ...
);