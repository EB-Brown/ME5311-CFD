function field = runge_kutta(field, c1, time_slope1, c2, time_slope2, dt)
    field = field + dt * (c1 * time_slope1 + c2 * time_slope2);
end