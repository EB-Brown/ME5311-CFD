function solver = laplacian_solver(k_mod, field)
    solver = idct2(k_mod.* dct2(field));
end