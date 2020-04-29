function dtheta_dt = get_dtheta_dt(...
    theta, ut_flux, vt_flux, dx, dy, x_num, y_num, ghost ...
)

% Convection term
dtheta_dt = -1 * ( ...
    (ut_flux(2:x_num + 1, :) - ut_flux(1:x_num, :)) / dx...
    + (vt_flux(:, 2:y_num + 1) - vt_flux(:, 1:y_num)) / dy...
);

gp1 = ghost + 1;
theta_dx_slice = theta(:, gp1:ghost + y_num);
theta_dy_slice = theta(gp1:ghost + x_num, :);

% Disperssion term
dtheta_dt = dtheta_dt + ( ...
    ( ...
        theta_dx_slice(1:x_num, :)  ...
        - 2 * theta_dx_slice(2:x_num + 1, :)  ...
        + theta_dx_slice(3:x_num + 2, :) ...
    ) / dx ^ 2 ...
    + ( ...
        theta_dy_slice(:, 1:y_num)  ...
        - 2 * theta_dy_slice(:, 2:y_num + 1)  ...
        + theta_dy_slice(:, 3:y_num + 2) ...
    ) / dy ^ 2 ...
);

end
