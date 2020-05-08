function [u_vel, v_vel, theta, heat_flux, time] = simulation_temp_flip( ...
    x_len, y_len, x_num, y_num, ...
    flip_time, total_time, cfl_target, ...
    prandtl, rayleigh, ...
    number_of_plots ...
)

%{
Purpose:
The objective of this simulation is to begin with repeating the
`simple_simulation.m` procedure until t = flip_time. Then the temperature
boundaries instantaneously switch walls and continue until t = total_time.

PARAMETERS
x_len: the length of the domain in the x direction
y_len: the length of the domain in the y direction
x_num: Number of grid cells in the x direction
y_num: Number of grid cells in the y direction
total_time: the amount of time the simulation will run for
flip_time: the time when the temperature domain will switch
cfl_target: The target CFL for controlling the time step size
prandtl: The Prandtl number for the simulation
rayleigh: The Rayleigh number for the simulation (Requires adjusting
    `cfl_target` accordingly to assure stability)
number_of_plots: The number of plots you would like to generate and save. A
    minimum of 1 is required to cover the initial condition. The time domain is
    divided into `end_time / (number_of_plots - 1)` chunks. A plot is generated
    after each chunk has been simulated.
%}

%%%%%%%%%%%%%%%%% Inputs %%%%%%%%%%%%%%%%%
% Temperature boundaries
left_wall_temp = 1;
right_wall_temp = 0;

% Number of ghost points
ghost = 1;

%%%%%%%%%%%%%%%%% Setup %%%%%%%%%%%%%%%%%
% Calculate space steps
dx = x_len/x_num;
dy = y_len/y_num;

% Calculate matrix for poisson solution
k_mod = get_k_mod(x_num, y_num, dx, dy);

% Initialize velocity as u(x,y) = v(x,y) = 0
[u_vel, v_vel] = initialize_velocity( ...
    x_num, y_num, dx, dy, k_mod, ghost ...
);

%Initialize temperature as as constant value
theta = zeros(x_num + 2 * ghost, y_num + 2 * ghost) + 0.5;

% Set temperature boundary condition
theta = temperature_boundaries( ...
    theta, x_num, y_num, ghost, left_wall_temp, right_wall_temp ...
);

% Initialize heat_flux and time structs
heat_flux.left = [];
heat_flux.right = [];
time = [];

% Divide plots between simpile simulation and flipped
total_iterations = number_of_plots - 1;
percent_simple = flip_time / total_time;
simple_iterations = floor(total_iterations * percent_simple);

%%%%%%%%%%%%%%%%% Output directory %%%%%%%%%%%%%%%%%
domain = "/x_len_" + num2str(x_len, "%.2g") ...
    + "__y_len_" + num2str(y_len, "%.2g");
domain = replace(replace(domain, ".", "p"), "+", "");

scaled_parms = "/pr_" + num2str(prandtl, "%.2e") ...
    + "__ra_" + num2str(rayleigh, "%.2e");
scaled_parms = replace(replace(scaled_parms, ".", "p"), "+", "");

output_dir = "plots/" + domain + scaled_parms + "/flip_temp_conditions";

%%%%%%%%%%%%%%%%% Simple Simulation %%%%%%%%%%%%%%%%%
[u_vel, v_vel, theta, heat_flux, time] = simulate_and_plot( ...
    x_len, y_len, x_num, y_num, dx, dy, k_mod, ghost, ...
    u_vel, v_vel, theta, ...
    left_wall_temp, right_wall_temp, heat_flux, ...
    cfl_target, time, 0, flip_time, ...
    simple_iterations, 0, ...
    prandtl, rayleigh, ...
    output_dir, 0 ...
);

%%%%%%%%%%%%%%%%% Flipped Temperature %%%%%%%%%%%%%%%%%
[u_vel, v_vel, theta, heat_flux, time] = simulate_and_plot( ...
    x_len, y_len, x_num, y_num, dx, dy, k_mod, ghost, ...
    u_vel, v_vel, theta, ...
    right_wall_temp, left_wall_temp, heat_flux, ... flipped temperature
    cfl_target, time, max(time), total_time, ...
    number_of_plots, simple_iterations + 1, ...
    prandtl, rayleigh, ...
    output_dir ...
);

end
