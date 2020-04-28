function [u_vel, v_vel] = velocity_boundaries(x_num, y_num, ghost, u_vel, v_vel)

    u_vel(ghost,:)=0; % Left - no penetration
    u_vel(x_num+ghost,:)=0; % right - no penetration

    v_vel(:,ghost)=0; % Bottom wall - no penetration
    v_vel(:,y_num+ghost)=0; % Top wall - no penetration

end