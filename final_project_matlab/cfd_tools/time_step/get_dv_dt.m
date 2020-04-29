function dv_dt = get_dv_dt( ...
    v_vel, theta, ...
    uv_flux, vv_flux, ...
    prandtl, rayleigh, ...
    dx, dy, x_num, y_num, ghost ...
)

% M-1 X N-2
% Convection Term
dv_dt = -1 * ( ...
    (uv_flux(2:x_num + 1, :) - uv_flux(1:x_num, :)) / dx...
    + (vv_flux(:, 2:y_num) - vv_flux(:, 1:y_num - 1)) / dy...
);

gp1 = ghost + 1;
right_idx = ghost + x_num;
roof_idx = ghost + y_num;
v_dx_slice = v_vel(ghost:right_idx + 1, gp1:roof_idx - 1);
v_dy_slice = v_vel(gp1:right_idx, ghost:roof_idx);

% Viscous term
dv_dt = dv_dt + prandtl * ( ...
    ( ...
        v_dx_slice(1:x_num, :)  ...
        - 2 * v_dx_slice(2:x_num + 1, :)  ...
        + v_dx_slice(3:x_num + 2, :) ...
    ) / dx ^ 2 ...
    + ( ...
        v_dy_slice(:, 1:y_num - 1)  ...
        - 2 * v_dy_slice(:, 2:y_num)  ...
        + v_dy_slice(:, 3:y_num + 1) ...
    ) / dy ^ 2 ...
);

% Bouyancy term
theta_slice = theta(gp1:right_idx, gp1:roof_idx)
dv_dt = dv_dt + (prandtl * rayleigh / 2) * ( ...
    theta(:, 2:y_num) + theta(:, 1:y_num - 1)  ...
);

end