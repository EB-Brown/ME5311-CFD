function [u_vel, v_vel] = time_march( ...
    u_vel, v_vel, dx, dy, dt, x_num, y_num, k_mod, ghost ...
)

gp1 = ghost + 1;
u_right_wall = x_num + ghost - 1;
v_right_wall = x_num + ghost - 2;
u_top_wall = y_num + ghost - 2;
v_top_wall = y_num + ghost - 1;

%%%%%%%%%%%%%%%% Initial increment %%%%%%%%%%%%%%%%
[du_dt, dv_dt] = get_time_slope( ...
    u_vel, v_vel, dx, dy, x_num, y_num, ghost ...
);

u_na = zeros(x_num + 2 * ghost, y_num + 2 * ghost);
u_na(gp1:u_right_wall - 1, ghost:u_top_wall) = runge_kutta( ...
    u_vel(gp1:u_right_wall - 1, ghost:u_top_wall), 8/15, du_dt, 0, 0, dt ...
);

v_na = zeros(x_num + 2 * ghost, y_num + 2 * ghost);
v_na(ghost:v_right_wall, gp1:v_top_wall - 1) = runge_kutta( ...
    v_vel(ghost:v_right_wall, gp1:v_top_wall - 1), 8/15, dv_dt, 0, 0, dt ...
);

[u_na, v_na] = remove_divergence( ...
    u_na, v_na, dx, dy, x_num, y_num, k_mod, ghost ...
);

%%%%%%%%%%%%%%%% Second increment %%%%%%%%%%%%%%%%
[du_na_dt, dv_na_dt] = get_time_slope( ...
    u_na, v_na, dx, dy, x_num, y_num, ghost ...
);

u_na1 = zeros(x_num + 2 * ghost, y_num + 2 * ghost);
u_na1(gp1:u_right_wall - 1, ghost:u_top_wall) = runge_kutta( ...
    u_na(gp1:u_right_wall - 1, ghost:u_top_wall), ...
    -17/60, ...
    du_dt, ...
    5/12, ...
    du_na_dt, ...
    dt ...
);

v_na1 = zeros(x_num + 2 * ghost, y_num + 2 * ghost);
v_na1(ghost:v_right_wall, gp1:v_top_wall - 1) = runge_kutta( ...
    v_na(ghost:v_right_wall, gp1:v_top_wall - 1), ...
    -17/60, ...
    dv_dt, ...
    5/12, ...
    dv_na_dt, ...
    dt...
);

[u_na1, v_na1] = remove_divergence( ...
    u_na1, v_na1, dx, dy, x_num, y_num, k_mod, ghost ...
);

%%%%%%%%%%%%%%%% Second increment %%%%%%%%%%%%%%%%
[du_na1_dt, dv_na1_dt] = get_time_slope( ...
    u_na1, v_na1, dx, dy, x_num, y_num, ghost ...
);

u_vel(gp1:u_right_wall - 1, ghost:u_top_wall) = runge_kutta( ...
    u_na1(gp1:u_right_wall - 1, ghost:u_top_wall), ...
    -5/12, ...
    du_na_dt, ...
    3/4, ...
    du_na1_dt, ...
    dt ...
);

v_vel(ghost:v_right_wall, gp1:v_top_wall - 1) = runge_kutta( ...
    v_na1(ghost:v_right_wall, gp1:v_top_wall - 1), ...
    -5/12, ...
    dv_na_dt, ...
    3/4, ...
    dv_na1_dt, ...
    dt...
);

[u_vel, v_vel] = remove_divergence( ...
    u_vel, v_vel, dx, dy, x_num, y_num, k_mod, ghost ...
);

end
