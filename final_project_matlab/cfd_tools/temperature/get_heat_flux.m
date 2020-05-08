function [left_flux, right_flux] = get_heat_flux( ...
    temperature, left_temp, right_temp, x_num, y_num, ghost ...
)

gp1 = ghost + 1;
theta = temperature(gp1:ghost+x_num, gp1:ghost + y_num);

left_flux = sum(temperature(1,:) - left_temp) / y_num;

along_right = temperature(ghost,gp1:ghost + y_num);
right_flux = sum(right_temp - temperature(x_num,:)) / y_num;

end
