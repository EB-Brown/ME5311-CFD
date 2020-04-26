function k_mod = get_k_mod(x_num, y_num, dx, dy)

    k_mod = zeros(x_num,y_num);

    for m=0:x_num-1
        for n=0:y_num-1

            k_mod(m + 1, n + 1) = 1 / ( ...
                (2 * cos(pi * m / x_num) - 2) / dx ^ 2 ...
                + (2 * cos(pi * n / y_num) - 2) / dy ^ 2 ...
            );

        end
    end

    k_mod(1, 1) = 0;

end