function [du_dt, dv_dt, dtheta_dt] = get_time_slope( ...
    u_vel, v_vel, theta, ...
    dx, dy, ...
    prandtl, rayleigh, ...
    x_num, y_num, ghost ...
)

[uu_flux, vu_flux, uv_flux, vv_flux, ut_flux, vt_flux] = get_flux( ...
    u_vel, v_vel, theta, x_num, y_num, ghost ...
);

% M-2 X N-1
du_dt = get_du_dt( ...
    u_vel, uu_flux, vu_flux, prandtl, dx, dy, x_num, y_num, ghost ...
);

% M-1 X N-2
dv_dt = get_dv_dt( ...
    v_vel, theta, ...
    uv_flux, vv_flux, ...
    prandtl, rayleigh, ...
    dx, dy, x_num, y_num, ghost ...
);

% M-1 X N-1
dtheta_dt = get_dtheta_dt(theta, ut_flux, vt_flux, dx, dy, x_num, y_num, ghost);

end