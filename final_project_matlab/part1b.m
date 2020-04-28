clear all
clc

% Number of points
x_num = 64;
y_num = 256;

x_len = 0.25;
y_len = 1;

ghost = 1; % number  of ghost cells

dx = x_len/x_num;
dy = y_len/y_num;

k_mod = get_k_mod(x_num, y_num, dx, dy);

% Fill u and v with random numbers
u_velocity = rand(x_num + 2 * ghost, y_num + 2 * ghost);
v_velocity = rand(x_num + 2 * ghost, y_num + 2 * ghost);

% Set boundary condition
[u_velocity, v_velocity] = velocity_boundaries( ...
    u_velocity, v_velocity, x_num, y_num, ghost ...
);

prior_divergence = velocity_divergence( ...
    u_velocity, v_velocity, dx, dy, x_num, y_num, ghost ...
);

max(max(prior_divergence))

[u_velocity, v_velocity] = remove_divergence( ...
    u_velocity, v_velocity, dx, dy, x_num, y_num, k_mod, ghost ...
);


after_divergence = velocity_divergence( ...
    u_velocity, v_velocity, dx, dy, x_num, y_num, ghost ...
);

max(max(after_divergence))
