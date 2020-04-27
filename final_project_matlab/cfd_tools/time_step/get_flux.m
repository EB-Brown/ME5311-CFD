function [uu_flux, vu_flux, uv_flux, vv_flux] = get_flux( ...
    u_vel, v_vel, x_num, y_num, ghost ...
)
% temperature
% ut_flux, vt_flux

% Trim off the ghost points
u = u_vel(ghost:x_num + ghost, ghost + 1:y_num + ghost); % M X N-1
v = v_vel(ghost + 1:x_num + ghost, ghost:y_num + ghost); % M-1 X N

% x_num: M-1
% y_num: N-1

u_right = x_num + 1;
v_right = x_num;

u_top = y_num; % Index of last row in U domain
v_top = y_num + 1;

uy_neg_2 = u_top - 1; % Index of second to last row in the U domain
vx_neg_2 = v_right - 1; % Index of second to last column in the V domain

% M-1 X N-1 -> dx -> M-2 X N-1
ux_shift = (u(2:u_right, :) + u(1:x_num, :)) / 2;

% Multi used for vu and uv flux. Kept as M X N and indexed on appropriate index
vx_shift = zeros(u_right, v_top);
vx_shift(2:x_num, :) = (v(2:v_right, :) + v(1:vx_neg_2, :)) / 2;

% Multi used for vu and uv flux. Kept as M X N and indexed on appropriate index
uy_shift = zeros(u_right, v_top);
uy_shift(:, 2:y_num) = (u(:, 2:u_top) + u(:, 1:uy_neg_2)) / 2;

vy_shift = (v(:, 2:v_top) + v(:, 1:y_num)) / 2;

uv = uy_shift.* vx_shift; % M X N

% U Flux
uu_flux = ux_shift.^ 2; % M-1 X N-1 -> dx -> M-2 X N-1
vu_flux = uv(2:x_num, :); % M-2 X N -> dy -> M-2 X N-1

% V Flux
uv_flux = uv(:, 2:y_num); % M X N-2 -> dx -> M-1 X N-2
vv_flux = vy_shift.^ 2; % M-1 X N-1 -> dy -> M-1 X N-2

end