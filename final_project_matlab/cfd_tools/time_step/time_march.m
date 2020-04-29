function [u_vel, v_vel, theta] = time_march( ...
    u_vel, v_vel, theta, ...
    dx, dy, dt, ...
    prandtl, rayleigh, ...
    x_num, y_num, ...
    k_mod, ghost ...
)

[vel_x_num, vel_y_num] = size(u_vel);

u_bott = ghost + 1;
v_bott = ghost;
vy_pls_1 = v_bott + 1;

u_left = ghost;
v_left = ghost + 1;
ux_pls_1 = u_left + 1;

right = x_num + ghost;
x_neg_2 = right - 1;

top = y_num + ghost;
y_neg_2 = top - 1;

%%%%%%%%%%%%%%%% Initial increment %%%%%%%%%%%%%%%%
[du_dt, dv_dt, dtheta_dt] = get_time_slope( ...
    u_vel, v_vel, theta, ...
    dx, dy, ...
    prandtl, rayleigh, ...
    x_num, y_num, ghost ...
);

theta_na = runge_kutta(theta, 8/15, dtheta_dt, 0, 0, dt);

u_na = zeros(vel_x_num, vel_y_num);
u_na(ux_pls_1:x_neg_2, u_bott:top) = runge_kutta( ...
    u_vel(ux_pls_1:x_neg_2, u_bott:top), 8/15, du_dt, 0, 0, dt ...
);

v_na = zeros(vel_x_num, vel_y_num);
v_na(v_left:right, vy_pls_1:y_neg_2) = runge_kutta( ...
    v_vel(v_left:right, vy_pls_1:y_neg_2), 8/15, dv_dt, 0, 0, dt ...
);

[u_na, v_na] = remove_divergence( ...
    u_na, v_na, dx, dy, x_num, y_num, k_mod, ghost ...
);

%%%%%%%%%%%%%%%% Second increment %%%%%%%%%%%%%%%%
[du_na_dt, dv_na_dt, dtheta_na_dt] = get_time_slope( ...
    u_na, v_na, theta_na, ...
    dx, dy, ...
    prandtl, rayleigh, ...
    x_num, y_num, ghost ...
);

theta_nb = runge_kutta(theta_na, -17/60, dtheta_dt, 5/12, dtheta_na_dt, dt);

u_nb = zeros(vel_x_num, vel_y_num);
u_nb(ux_pls_1:x_neg_2, u_bott:top) = runge_kutta( ...
    u_na(ux_pls_1:x_neg_2, u_bott:top), -17/60, du_dt, 5/12, du_na_dt, dt ...
);

v_nb = zeros(vel_x_num, vel_y_num);
v_nb(v_left:right, vy_pls_1:y_neg_2) = runge_kutta( ...
    v_na(v_left:right, vy_pls_1:y_neg_2), -17/60, dv_dt, 5/12, dv_na_dt, dt ...
);

[u_nb, v_nb] = remove_divergence( ...
    u_nb, v_nb, dx, dy, x_num, y_num, k_mod, ghost ...
);

%%%%%%%%%%%%%%%% Second increment %%%%%%%%%%%%%%%%
[du_nb_dt, dv_nb_dt, dtheta_nb_dt] = get_time_slope( ...
    u_nb, v_nb, theta_nb, ...
    dx, dy, ...
    prandtl, rayleigh, ...
    x_num, y_num, ghost ...
);

theta = runge_kutta(theta_nb, -5/12, dtheta_na_dt, 3/4, dtheta_nb_dt, dt);

u_vel(ux_pls_1:x_neg_2, u_bott:top) = runge_kutta( ...
    u_nb(ux_pls_1:x_neg_2, u_bott:top), -5/12, du_na_dt, 3/4, du_nb_dt, dt ...
);

v_vel(v_left:right, vy_pls_1:y_neg_2) = runge_kutta( ...
    v_nb(v_left:right, vy_pls_1:y_neg_2), -5/12, dv_na_dt, 3/4, dv_nb_dt, dt...
);

[u_vel, v_vel] = remove_divergence( ...
    u_vel, v_vel, dx, dy, x_num, y_num, k_mod, ghost ...
);

end
