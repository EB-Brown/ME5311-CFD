clear all
clc

%{
This script bulk processes the three simulation types listed below:
    simulation_simple.m - Standard boundary conditions
    simulation_sine_temp.m - Temperature initialized as a sine wave
    simulation_temp_flip.m - Temperature boundaries flip mid-simulation

Each of the above functions will generate plots for the fluid profile and the
heat flus along the walls.

For more information regarding these simulation structures, please check the
documentation in each script.

This batch script will create unique combinations of the domains, prandtl
numbers, and rayleigh numbers, then run simulations on each combination storing
the results under the "plots" folder under the following directory structure:

    domain/prandtl_rayleigh/simulation
%}

%%%%%%%%%%% Inputs %%%%%%%%%%%%%
domains = [ ... [x length, y length, x grid cells, y grid cells]
    [0.75, 2, 75, 200]; ...
    [2, 2, 150, 150]; ...
    [2, 1, 200, 100]; ...
];

prandtl_numbers = [0.7, 4];
rayleigh_numbers = [1e4, 1e8];

end_time = 0.01; % Total run time for simple and sine, flip time for temp_flip
cfl_target = 0.3;
number_of_plots = 10;

%%%%%%%%%%% Setup %%%%%%%%%%%
simulations = zeros( ...
    6, length(domains) * length(prandtl_numbers) * length(rayleigh_numbers) ...
);

m = 0;
for domain=domains'
    for prandtl=prandtl_numbers
        for rayleigh=rayleigh_numbers
            m = m + 1;
            simulations(:, m) = [domain; prandtl; rayleigh];
        end
    end
end

total_sims = size(simulations, 2);

flip_time = end_time;
total_time_for_flip = 3 * end_time;
total_plots_for_flip = 3 * number_of_plots;

%%%%%%%%%%% Running simulations %%%%%%%%%%%
progress = 0;
for simulation=simulations

    x_len = simulation(1);
    y_len = simulation(2);
    x_num = simulation(3);
    y_num = simulation(4);
    prandtl = simulation(5);
    rayleigh = simulation(6);

    simulation_simple( ...
        x_len, y_len, x_num, y_num, ...
        end_time, cfl_target, ...
        prandtl, rayleigh, ...
        number_of_plots ...
    );

    simulation_sine_temp( ...
        x_len, y_len, x_num, y_num, ...
        end_time, cfl_target, ...
        prandtl, rayleigh, ...
        number_of_plots ...
    );

    simulation_temp_flip( ...
        x_len, y_len, x_num, y_num, ...
        flip_time, total_time_for_flip, cfl_target, ...
        prandtl, rayleigh, ...
        total_plots_for_flip ...
    );

    clc
    progress = progress + 1;
    fprintf("Progress : %.2d \n", 100 * progress/total_sims)

end
