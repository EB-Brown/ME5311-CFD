function field = pressure_dy(field, dy, x_num, y_num, ghost)

field = ( ...
    field(ghost + 1:x_num + ghost, ghost + 2:y_num + ghost) ...
    - field(ghost + 1:x_num + ghost, ghost + 1:y_num + ghost - 1) ...
) / dy;

end