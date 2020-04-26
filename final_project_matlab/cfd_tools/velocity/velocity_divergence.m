function vel_div = velocity_divergence( ...
    u_vel, v_vel, dx, dy, x_num, y_num, ghost ...
)

part_x = velocity_dx(u_vel, dx, x_num, y_num, ghost);
part_y = velocity_dy(v_vel, dy, x_num, y_num, ghost);
vel_div = part_x + part_y;

end
