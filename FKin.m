function [T] = FKin(dhparams, theta)
[dof,~] = size(theta);
T = eye(4);
for i = 1:dof
    a  = dhparams(i,1);
    al = dhparams(i,2);
    if dhparams(i,4) == 1
        d = dhparams(i,3);
        th = theta(i);
    else
        th = dhparams(i,3);
        d  = theta(i);
    end
%     Ti = [cos(th) -sin(th) 0 a;
%         sin(th)*cos(al) cos(th)*cos(al) -sin(al) -d*sin(al);
%         sin(th)*sin(al) cos(th)*sin(al) cos(al) d*cos(al);
%         0 0 0 1];
    Ti = [cos(th)  -sin(th)*cos(al) sin(th)*sin(al)  a*cos(th);
              sin(th)  cos(th)*cos(al)  -cos(th)*sin(al) a*sin(th);
              0        sin(al)          cos(al)          d;
              0                0                0        1];
    T = T * Ti;
end
end