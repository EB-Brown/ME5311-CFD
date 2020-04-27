clear all
clc

%%%%%%%%%%%%%%%%% Inputs %%%%%%%%%%%%%%%%%
% Number of points
x_num = 64;
y_num = 256;
cfl_target = 1.19; % Initial time step calculated from CFL
time_iterations = 100; % Minimum number of time_iterations

% Number of simulations to run with decreasing time steps
simulation_iterations = 10;

% Controls tested dt. ex: .5 will make final tested dt = 0.5 initial dt
final_dt_scale = 0.5;

x_len = 0.25;
y_len = 1;

ghost = 1; % number  of ghost cells

%%%%%%%%%%%%%%%%% Setup %%%%%%%%%%%%%%%%%

% find boundaries
top_wall = ghost + y_num;
right_wall = ghost + x_num;

% Calculate space steps
dx = x_len/x_num;
dy = y_len/y_num;

% Calculate matrix for poisson solution
k_mod = get_k_mod(x_num, y_num, dx, dy);

% Fill u and v with random numbers
u_velocity = rand(x_num + 2 * ghost, y_num + 2 * ghost);
v_velocity = rand(x_num + 2 * ghost, y_num + 2 * ghost);

% Set boundary condition
[u_velocity, v_velocity] = boundaries( ...
    x_num, y_num, ghost, u_velocity, v_velocity ...
);

% Remove Divergence
[u_velocity, v_velocity] = remove_divergence( ...
    u_velocity, v_velocity, dx, dy, x_num, y_num, k_mod, ghost ...
);

% Calculate dt size
initial_dt = cfl_target / ( ...
    max(max(abs(u_velocity)))/dx + max(max(abs(v_velocity)))/dy ...
);

end_time = time_iterations * initial_dt; % Seconds

% Calculate initial kinetic energy
initial_kinetic_energy = get_kinetic_energy( ...
    u_velocity, v_velocity, top_wall, right_wall, ghost ...
);
initial_ke_sum = sum(sum(initial_kinetic_energy));

% An array for tracking kinetic energy with respect to dt
kinetic_trend = zeros(3, simulation_iterations);

%%%%%%%%%%%%%%%%% Time Integration %%%%%%%%%%%%%%%%%

for n=1:simulation_iterations

    factor = (1 - final_dt_scale * (n - 1)/(simulation_iterations - 1));
    dt = initial_dt * factor;

    % Run simulation
    [u_vel, v_vel] = simulate( ...
        u_velocity, v_velocity, dx, dy, ...
        dt, 0, end_time, ...
        x_num, y_num, ghost, ...
        k_mod  ...
    );

    % Calculate kinetic energy
    kinetic_energy = get_kinetic_energy( ...
        u_vel, v_vel, top_wall, right_wall, ghost ...
    );
    ke_sum = sum(sum(kinetic_energy));

    kinetic_trend(1, n) = dt;
    kinetic_trend(2, n) = abs(ke_sum - initial_ke_sum);
    kinetic_trend(3, n) = ke_sum;

    fprintf("Progress: %.2d\n", 100 * n / simulation_iterations);

end

fig = kinetic_energy_plot( ...
    kinetic_trend, initial_ke_sum, "kinetic_energy_convergence.png" ...
);

disp("Simulations Complete")
