function [u_vel, v_vel, theta, heat_flux, time] = simulation_sine_temp( ...
    x_len, y_len, x_num, y_num, ...
    end_time, cfl_target, ...
    prandtl, rayleigh, ...
    number_of_plots ...
)

%{
Procedure:
This simulation initializes the temperature profile with a quarter sin wave. The
left wall has a fixed temperature of 1 and the right wall has a fixed temperatre
of 0. The sin wave is defined as sin(pi* x_array / (2 * x_length)) so that the
fluid temperature toward the left wall approaches 0 and the fluid temperature
toward the right wall approaches 1.

PARAMETERS
x_len: the length of the domain in the x direction
y_len: the length of the domain in the y direction
x_num: Number of grid cells in the x direction
y_num: Number of grid cells in the y direction
end_time: the amount of time the simulation will run for
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

% Initial horizontal temperature profile
top_wall = ghost + y_num;
theta_x = dx/2:dx:x_len-dx/2;
init_temp = sin(pi*theta_x / (2*x_len));

% Assign profile to each row of temperature
gp1 = ghost + 1;
right_wall = ghost + x_num;
theta = zeros(x_num + 2 * ghost, y_num + 2 * ghost);
for row=ghost+1:top_wall
    theta(gp1:right_wall, row) = init_temp;
end

% Set temperature boundary condition
theta = temperature_boundaries( ...
    theta, x_num, y_num, ghost, left_wall_temp, right_wall_temp ...
);

% Initialize heat_flux and time structs
heat_flux.left = [];
heat_flux.right = [];
time = [];

%%%%%%%%%%%%%%%%% Output directory %%%%%%%%%%%%%%%%%
domain = "/x_len_" + num2str(x_len, "%.2g") ...
    + "__y_len_" + num2str(y_len, "%.2g");
domain = replace(replace(domain, ".", "p"), "+", "");

scaled_parms = "/pr_" + num2str(prandtl, "%.2e") ...
    + "__ra_" + num2str(rayleigh, "%.2e");
scaled_parms = replace(replace(scaled_parms, ".", "p"), "+", "");

output_dir = "plots/" + domain + scaled_parms + "/sine_temp_conditions";

%%%%%%%%%%%%%%%%% Simulation %%%%%%%%%%%%%%%%%
[u_vel, v_vel, theta, heat_flux, time] = simulate_and_plot( ...
    x_len, y_len, x_num, y_num, dx, dy, k_mod, ghost, ...
    u_vel, v_vel, theta, ...
    left_wall_temp, right_wall_temp, heat_flux, ...
    cfl_target, time, 0, end_time, ...
    number_of_plots, 0, ...
    prandtl, rayleigh, ...
    output_dir ...
);

end
