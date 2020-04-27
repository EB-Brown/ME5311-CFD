function kinetic_energy = get_kinetic_energy( ...
    u_vel, v_vel, top_wall, right_wall, ghost ...
)
u = u_vel(ghost:right_wall, ghost:top_wall);
v = v_vel(ghost:right_wall, ghost:top_wall);

kinetic_energy = sum(sum(u.^2 + v.^2));

end