clear all
clc

%%%%%%%%%%%%%%%%% Inputs %%%%%%%%%%%%%%%%%
% Number of points
x_num = 64;
y_num = 256;
cfl_target = 1.1;
time_iterations = 1000;

x_len = 0.25;
y_len = 1;

ghost = 1; % number  of ghost cells

%%%%%%%%%%%%%%%%% Setup %%%%%%%%%%%%%%%%%

top_wall = ghost + y_num;
right_wall = ghost + x_num;

dx = x_len/x_num;
dy = y_len/y_num;

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

%%%%%%%%%%%%%%%%% Time Integration %%%%%%%%%%%%%%%%%
time = 0;
subplot(2,2,1)
contour(u_velocity,100)
colorbar()
title("U Velocity")
subplot(2,2,3)
contour(v_velocity,100)
title("V Velocity")
colorbar()
% pause(0.01)

dt = cfl_target / ( ...
    max(max(abs(u_velocity)))/dx + max(max(abs(v_velocity)))/dy ...
);
initial_kinetic_energy = get_kinetic_energy( ...
    u_velocity, v_velocity, top_wall, right_wall, ghost ...
);
kinetic_trend = zeros(2,time_iterations + 1);
kinetic_trend(2, 1) = sum(sum(initial_kinetic_energy));

for n=1:time_iterations
    time = time + dt;

    [u_velocity, v_velocity] = time_march( ...
        u_velocity, v_velocity, dx, dy, dt, x_num, y_num, k_mod, ghost ...
    );

    dt = get_next_dt(cfl_target, u_velocity, v_velocity, dx, dy, dt);

    kinetic_energy = get_kinetic_energy( ...
        u_velocity, v_velocity, top_wall, right_wall, ghost ...
    );
    kinetic_trend(1, n + 1) = time;
    kinetic_trend(2, n + 1) = sum(sum(kinetic_energy));
    % kinetic_change = kinetic_energy - initial_kinetic_energy;

    %{
    diverg = velocity_divergence( ...
        u_velocity, v_velocity, dx, dy, x_num, y_num, ghost ...
    );

    fprintf("iteration = %g \n", n);
    fprintf("time step = %g \n", dt);
    fprintf("time = %g \n", time);
    fprintf("max divergence = %.2e \n", max(max(abs(diverg))));
    fprintf("max U = %.2f \n", max(max(abs(u_velocity))));
    fprintf("max V = %.2f \n", max(max(abs(v_velocity))));
    fprintf("Kinetic Change = %.2f \n", kinetic_change);
    fprintf("\n")
    %}
end

subplot(2,3,1)
contour(u_velocity(ghost:ghost + x_num, ghost + 1: ghost + y_num),100)
colorbar()
title("U Velocity")
caxis([-1, 1])
subplot(2,3,4)
contour(v_velocity(ghost + 1:ghost + x_num, ghost: ghost + y_num),100)
title("V Velocity")
colorbar()
caxis([-1, 1])
subplot(2,3,2)
plot(kinetic_trend(1, 1:n+1), kinetic_trend(2, 1:n+1));
title("Total Kinetic Energy")
subplot(2,3,5)
plot(kinetic_trend(1, 1:n), diff(kinetic_trend(2, 1:n+1)));
title("Change in Kinetic Energy")
subplot(2,3,3)
contour(kinetic_energy,100)
title("Kinetic Energy")
colorbar()

disp("Done")

% pause(0.001)

%end
