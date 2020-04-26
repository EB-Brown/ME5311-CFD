function field = velocity_dx(field, dx, x_num, y_num, ghost)

field = ( ...
    field(ghost + 1:x_num + ghost, ghost + 1:y_num + ghost) ...
    - field(ghost:x_num + ghost - 1, ghost + 1:y_num + ghost) ...
) / dx;

end