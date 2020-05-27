function [T, ll, J] = claculate_dhparam_and_jacobian()
% the parameters & link limits for 6 dof robot as described in SVD paper
    dhparam = [0 pi/2 537.5 1;32.5 -pi/2 0 1;0 0 0 0;43.8 pi/2 73.6 1;0 -pi/2 0 1;0 0 208.2 1];
    link_limits(1,1:2) = [-pi/4 pi/4];
    link_limits(2,1:2) = [-pi/4 pi/4];
    link_limits(3,1:2) = [0 0.5];
    link_limits(4,1:2) = [-pi/2 pi/2];
    link_limits(5,1:2) = [-pi/2 pi/2];
    link_limits(6,1:2) = [-pi/2 pi/2];
    disp = [pi/2 pi/2 956.6 0 0 0]';
    disp(:,2) = disp;
    ll = link_limits + disp;
    
% calculate the associated jacobian of the spatial twist
    syms th1 th2 th3 th4 th5 th6;
    assume([th1 th2 th3 th4 th5 th6], 'real');
    
%     S = [1 0 0 0 537.5 0;...
%         0 -1 0 570 0 0;...
%         0 0 0 0 -1 0;...
%         1 0 0 0 (570+43.8) (956.6+th3+73.6);...
%         0 -1 0 (570+43.8) 0 0;...
%         0 -1 0 (570+43.8) 0 0;]';
    
    joint_vars = [th1 th2 th3 th4 th5 th6];
    
    T = FKin(dhparam, joint_vars');
    r1 = [diff(T(1,4),th1) diff(T(1,4),th2) diff(T(1,4),th3) diff(T(1,4),th4) diff(T(1,4),th5) diff(T(1,4),th6)];
    r2 = [diff(T(2,4),th1) diff(T(2,4),th2) diff(T(2,4),th3) diff(T(2,4),th4) diff(T(2,4),th5) diff(T(2,4),th6)];
    r3 = [diff(T(3,4),th1) diff(T(3,4),th2) diff(T(3,4),th3) diff(T(3,4),th4) diff(T(3,4),th5) diff(T(3,4),th6)];
    J = [r1;r2;r3];
  
end

function [Tinv] = Transformation_Matrix_Inverse(T)
    p = T(1:3,4);
    R = T(1:3,1:3);
    R_inverse = R';
    p_inverse = -1 * R_inverse * p;
    Tinv = [R_inverse p_inverse;zeros(1,3) 1];
end