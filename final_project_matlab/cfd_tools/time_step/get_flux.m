function [uu_flux, vu_flux, uv_flux, vv_flux] = get_flux( ...
    u_vel, v_vel, x_num, y_num, ghost ...
)
% temperature
% ut_flux, vt_flux

% Trim off the ghost points
u = u_vel(ghost:x_num + ghost - 1, ghost:y_num + ghost - 2); % M X N-1
v = v_vel(ghost:x_num + ghost - 2, ghost:y_num + ghost - 1); % M-1 X N

ux_neg_2 = x_num - 1; % Index of second to last column in the U domain
vx_neg_1 = ux_neg_2; % Index of last column in V domain
vx_neg_2 = vx_neg_1 - 1; % Index of second to last column in the V domain

vy_neg_2 = y_num - 1; % Index of second to last row in the V domain
uy_neg_1 = vy_neg_2; % Index of last row in U domain
uy_neg_2 = uy_neg_1 - 1; % Index of second to last row in the U domain

% M-1 X N-1 -> dx -> M-2 X N-1
ux_shift = (u(2:x_num, :) + u(1:ux_neg_2, :)) / 2;

% M-2 X N -> dy -> M-2 X N-1
vx_shift = zeros(x_num, y_num);
vx_shift(2:vx_neg_1, :) = (v(2:vx_neg_1, :) + v(1:vx_neg_2, :)) / 2;

% M-2 X N -> dy -> M-2 X N-1
uy_shift = zeros(x_num, y_num);
uy_shift(:, 2:uy_neg_1) = (u(:, 2:uy_neg_1) + u(:, 1:uy_neg_2)) / 2;

vy_shift = (v(:, 2:y_num) + v(:, 1:vy_neg_2)) / 2;

uv = uy_shift.* vx_shift; % M X N

% U Flux
uu_flux = ux_shift.^ 2; % M-1 X N-1 -> dx -> M-2 X N-1
vu_flux = uv(2:ux_neg_2, :); % M-2 X N -> dy -> M-2 X N-1

% V Flux
uv_flux = uv(:, 2:vy_neg_2); % M X N-2 -> dx -> M-1 X N-2
vv_flux = vy_shift.^ 2; % M-1 X N-1 -> dy -> M-1 X N-2

end