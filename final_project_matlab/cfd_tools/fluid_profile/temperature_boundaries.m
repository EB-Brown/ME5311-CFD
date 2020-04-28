function temperature = temperature_boundaries( ...
    temperature, x_num, y_num, ghost, left_condition ...
)

% Averages are manipulated due to staggered grid

% Top and bottom wall partial temp / partial y = 0
temperature(:, ghost) = temperature(:, ghost + 1); % Bottom
temperature(:, ghost + y_num + 1) = temperature(:, ghost + y_num); % Top

% Use the given left wall condition
temperature(ghost, :) = 2 * left_condition - temperature(ghost + 1, :);

% Right wall = 0
temperature(ghost + x_num + 1, :) = -1 * temperature(ghost + x_num, :);

end

