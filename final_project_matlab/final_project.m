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
subplot(2,1,1)
contour(u_velocity,100)
colorbar()
subplot(2,1,2)
contour(v_velocity,100)
colorbar()
pause(0.01)

dt = cfl_target / ( ...
    max(max(abs(u_velocity)))/dx + max(max(abs(v_velocity)))/dy ...
);

for n=1:time_iterations
    [u_velocity, v_velocity] = time_march( ...
        u_velocity, v_velocity, dx, dy, dt, x_num, y_num, k_mod, ghost ...
    );

    diverg = velocity_divergence( ...
        u_velocity, v_velocity, dx, dy, x_num, y_num, ghost ...
    );

    time = time + dt;

    fprintf("iteration = %g \n", n);
    fprintf("time step = %g \n", dt);
    fprintf("time = %g \n", time);
    fprintf("max divergence = %.2e \n", max(max(abs(diverg))));
    fprintf("max U = %.2f \n", max(max(abs(u_velocity))));
    fprintf("max V = %.2f \n", max(max(abs(v_velocity))));
    fprintf("\n")

    subplot(2,1,1)
    contour(u_velocity,100)
    colorbar()
    title("U Velocity")
    subplot(2,1,2)
    contour(v_velocity,100)
    title("V Velocity")
    colorbar()
    
    dt = get_next_dt(cfl_target, u_velocity, v_velocity, dx, dy, dt);

    pause(0.01)

end
