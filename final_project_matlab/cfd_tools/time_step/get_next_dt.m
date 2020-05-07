function dt = get_next_dt(cfl_target, u_vel, v_vel, prandtl, dx, dy, dt)

max_u = max(max(abs(u_vel))) ;
max_v = max(max(abs(v_vel))) ;
if (max_u == 0) && (max_v == 0)
        max_u = 1;
        max_v = 1;
end

cfl_convection = dt * (max_u / dx + max_v / dy) ;
cfl_diffusion = prandtl * dt * (1 / dx ^ 2 + 1 / dy ^ 2);
cfl_max = max(cfl_diffusion, cfl_convection);

dt = 0.5 * dt * (cfl_target / cfl_max + 1);

end
