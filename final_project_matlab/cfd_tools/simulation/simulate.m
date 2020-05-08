function [u_vel, v_vel, theta, heat_flux, time] = simulate( ...
        u_vel, v_vel, theta, heat_flux, dx, dy, ...
        dt, cfl, time, start_time, end_time, ...
        prandtl, rayleigh, ...
        x_num, y_num, left_temp, right_temp, ...
        ghost, k_mod  ...
)

if cfl % CFL is defined, use dynamic dt
    max_u = max(max(abs(u_vel)));
    max_v = max(max(abs(v_vel)));
    if (max_u == 0) && (max_v == 0)
        max_u = 1;
        max_v = 1;
    end

    dt_convection = cfl / (max_u/dx + max_v/dy);
    dt_diffusion = cfl / (prandtl * (1/dx^2 + 1/dy^2));
    dt = min(dt_convection, dt_diffusion);

end

current_time = start_time;
while true

    time = [time; current_time];

    [left_flux, right_flux] = get_heat_flux( ...
        theta, left_temp, right_temp, x_num, y_num, ghost ...
    );

    heat_flux.left = [heat_flux.left; left_flux];
    heat_flux.right = [heat_flux.right; right_flux];

    [u_vel, v_vel, theta] = time_march( ...
        u_vel, v_vel, theta, ...
        dx, dy, dt, ...
        prandtl, rayleigh, ...
        x_num, y_num, left_temp, right_temp,...
        k_mod, ghost ...
    );

    current_time = current_time + dt;
    if current_time >= end_time
        break
    end

    if cfl % CFL is defined, use dynamic dt
        dt = get_next_dt(cfl, u_vel, v_vel, prandtl, dx, dy, dt);
    end

end

end