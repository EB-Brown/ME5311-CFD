function dt = get_next_dt(cfl_target, u_vel, v_vel, dx, dy, dt)

cfl_max = dt * (max(max(abs(u_vel))) / dx + max(max(abs(v_vel))) / dy) ;
dt = 0.5 * dt * (cfl_target / cfl_max + 1);

end
