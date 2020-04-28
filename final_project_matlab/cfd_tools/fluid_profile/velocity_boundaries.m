function [u_vel, v_vel] = velocity_boundaries(u_vel, v_vel, x_num, y_num, ghost)

u_vel(ghost, :) = 0; % Left - no penetration
u_vel(x_num + ghost, :) = 0; % right - no penetration
u_vel(:, ghost) = -1 * u_vel(:, ghost + 1); % Bottom wall - average row 1 & 2= 0
u_vel(:, y_num + ghost + 1) = -1 * u_vel(:, y_num + ghost); % top average = 0

v_vel(:, ghost) = 0; % Bottom wall - no penetration
v_vel(:, y_num + ghost) = 0; % Top wall - no penetration
v_vel(ghost, :) = -1 * v_vel(ghost + 1, :); % right wall - average col 1 & 2= 0
v_vel(x_num + ghost + 1, :) = -1 * v_vel(x_num + ghost, :); % left average = 0

end