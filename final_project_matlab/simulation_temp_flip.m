clear all
clc

%{
Purpose:
The objective of this simulation is to begin with repeating the
`simple_simulation.m` procedure until t = 0.01. Then the temperature boundaries
instantaneously switch walls and continue until t = 0.03.

Hypothesis:
At t = 0.01, the flow will be steady in the clockwise direction. When the
temperature boundaries switch, the flow will move will be pulled in the opposite
direction. A vortex sheet will develop due to the oposing velocities and cause a
chaotic flow. The flow will settle with a the top part of the domain being hot
and the bottom portion being cold as minor natural convection exists on the walls
in the counter-clockwise direction.

Observations:

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
total_time = 0.03;
flip_time = 0.01;
left_wall_temperature = 1;
right_wall_temperature = 0;

% Dynamic time steps are calculated using the relaxation method
cfl_target = 0.6;

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

%Initialize temperature as as constant value
theta = zeros(x_num + 2 * ghost, y_num + 2 * ghost) + 0.5;

% Set temperature boundary condition
theta = temperature_boundaries( ...
    theta, ...
    x_num, ...
    y_num, ...
    ghost, ...
    left_wall_temperature, ...
    right_wall_temperature ...
);

heat_flux.left = [];
heat_flux.right = [];
time = [];

% Divide plots between simpile simulation and flipped
total_iterations = number_of_plots - 1;
percent_simple = flip_time / total_time;

simple_iterations = floor(total_iterations * percent_simple);
simple_time_chunks = flip_time / simple_iterations;

flip_iterations = total_iterations - simple_iterations;
flip_time_chunks = (total_time - flip_time) / flip_iterations;

%%%%%%%%%%%%%%%%% Simulation %%%%%%%%%%%%%%%%%

output_dir = "plots/flipped_temperature_conditions";
scaled_parms = "/"+ "pr_" + num2str(prandtl, "%.2e") ...
    + "__ra_" + num2str(rayleigh, "%.2e");
output_dir = output_dir + replace(replace(scaled_parms, ".", "_"), "+", "");

%%%%%%%%%%%%%%%%% Simple Simulation %%%%%%%%%%%%%%%%%
[u_velocity, v_velocity, theta, heat_flux, time] = simulate_and_plot( ...
    x_len, y_len, x_num, y_num, dx, dy, k_mod, ghost, ...
    u_velocity, v_velocity, theta, ...
    left_wall_temperature, right_wall_temperature, heat_flux, ...
    cfl_target, time, 0, flip_time, ...
    simple_iterations, 0, ...
    prandtl, rayleigh, ...
    output_dir, 0 ...
);

%%%%%%%%%%%%%%%%% Flipped Temperature %%%%%%%%%%%%%%%%%
[u_velocity, v_velocity, theta, heat_flux, time] = simulate_and_plot( ...
    x_len, y_len, x_num, y_num, dx, dy, k_mod, ghost, ...
    u_velocity, v_velocity, theta, ...
    right_wall_temperature, left_wall_temperature, heat_flux, ... flipped temp
    cfl_target, time, max(time), total_time, ...
    number_of_plots, simple_iterations, ...
    prandtl, rayleigh, ...
    output_dir ...
);

