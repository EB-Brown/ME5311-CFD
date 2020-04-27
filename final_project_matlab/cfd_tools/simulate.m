function [u_vel, v_vel, theta] = simulate( ...
        u_vel, v_vel, theta, dx, dy,  ...
        dt, cfl, end_time, ...
        x_num, y_num, ghost, ...
        k_mod  ...
)

if cfl % CFL is defined, use dynamic dt
    dt = cfl / (max(max(abs(u_vel)))/dx + max(max(abs(v_vel)))/dy);
end

time = 0;
while true

    [u_vel, v_vel, theta] = time_march( ...
        u_vel, v_vel, theta, dx, dy, dt, x_num, y_num, k_mod, ghost ...
    );

    time = time + dt;
    if time >= end_time
        break
    end

    if cfl % CFL is defined, use dynamic dt
        dt = get_next_dt(cfl_target, u_velocity, v_velocity, dx, dy, dt);
    end

end

end