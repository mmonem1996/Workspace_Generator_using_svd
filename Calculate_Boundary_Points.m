function [xyz] = Calculate_Boundary_Points(r, max_points, slice_curve_count, sl_curv_index)
    % the parameters & link limits for 6 dof robot as described in SVD paper
    [T, ll, J] = claculate_dhparam_and_jacobian();
    % obtaining initial values for the slice curves
    [xyz,Theta_min, Theta_max] = find_slice_curves(slice_curve_count, 3, 2, T, ll);
    % we will do the calculations via symbolic objects & inline functions
    p = matlabFunction(T(1:3,4));   
    % initial searching direction
    d1 = [-1 -0.5]';
    d = d1 / norm(d1);
    % getting the slice curve value
    curve_z = xyz(sl_curv_index,1);
    q1 = Theta_min(sl_curv_index,:)';
    qmax = Theta_max(sl_curv_index,:)';
    % q1 = qmax;
    q = q1;
    Jf = matlabFunction(J(1:2,1:5));
    p1 = p(q(1), q(2), q(3), q(4), q(5));
    % this is to visualize curve generation
    figure
    hold on
    plot(p1(1), p1(2),'or','MarkerSize',5,'MarkerFaceColor','r')
    pause(0.1)
    %
    fprintf('point 1 = (%0.2f, %0.2f, %0.2f)\n', p1(1),p1(2),p1(3));
  
    xyz = p1; %[p1 p2]
    p2 = p1;
    counter = 2;
    dist_from_p0 = 1.5*r;

    % prepare functions for optimization
    sing_func = get_sing_func(Jf, ll(1:5,1), ll(1:5,2));
    opt = optimset('display','off','MaxIter',3000,'MaxFunEvals',3000);
    
    w = 0;
    
    while (dist_from_p0 > 0.9 * r)
        estimation_objective_func = objective_func(r, d);
        nlcon = boundary_nolcon(r, Jf, q, p, 3,curve_z, ll(1:5,1), ll(1:5,2));
        w = fmincon(estimation_objective_func, w, [], [], [], [], [], [], nlcon, opt);
        % get the value of q associated with the optimum angle w
        v = r * [cos(w) sin(w)]';
        Ji = Jf(q(1), q(2), q(3), q(4), q(5));
        R = double(colspace(sym(Ji')));
        M = Ji * R;
        dq = R * (M \ v);
        q = q + dq;
        p2_prev = p2;
        p2 = p(q(1), q(2), q(3), q(4), q(5));
        % this is to visualize curve generation
        plot(p2(1), p2(2) ,'or','MarkerSize',5,'MarkerFaceColor','r')
        pause(0.1)
        %
        % q is the approximate coordinates of the boundary point
        % here we need to optimize q to find the accurate boundary point
        sing_nlcon = get_sing_nlcon(p, 3, curve_z, r, p2);
        q = fmincon(sing_func, q, [], [], [], [], ll(1:5,1), ll(1:5,2), sing_nlcon, opt);
        p2 = p(q(1), q(2), q(3), q(4), q(5));
        % ssval = sing_func(q)
        %%%%%%
        fprintf('point %d = (%0.2f, %0.2f, %0.2f)\n', counter, p(q(1), q(2), q(3), q(4), q(5))');
        % calculating new search direction d
        dpn = (p2 - p2_prev);
        % dpn = dpn / norm(dpn);
        d = dpn(1:2,1);
        d = d / norm(d)
        % calculate the distance between the currently calculated boundary  
        % point and the starting point, if the distance is small enough 
        % then the curve is closed 
        dist_from_p0 = norm(p2 - p1);
        counter = counter + 1;
        xyz = [xyz p2];
        if counter > max_points
            break
        end
    end
    xyz = xyz';
end

function [ob_func] = objective_func(r, d)
    function [f] = func(w)
        % the v vector is assumed to be a planar victor with constant
        % radius r, the test is whether the current value of the angle w
        % correspond to the maximum value of the dot product between v and
        % the search direction vector d
        v = r * [cos(w) sin(w)]';
        f = -1 * dot(d, v);
    end
    ob_func = @func;
end

function [nlcon] = boundary_nolcon(r, Jf, q, p, pi, p_con_val, ql, qu)%p, pi, p_con_val,
    function [c, ceq] = nl_func(w)
        v = r * [cos(w) sin(w)]';
        J = Jf(q(1), q(2), q(3), q(4), q(5));
        R = double(colspace(sym(J')));
        M = J * R;
        dq = R * (M \ v);
        A = (dq + q) - qu;
        B = ql - (dq + q);
%         pval = p(q(1), q(2), q(3), q(4), q(5));
%         C = pval(pi)-1.01*p_con_val;
%         D = 0.99 * p_con_val - pval(pi);
%         D = p_con_val - pval(pi);
        c = [A;B];
        ceq =[];% D;
    end
    nlcon = @nl_func;
end

function [sing_func] = get_sing_func(J, ql, qu)
    function [sing_val] = evaluate_singular_value(q)
        Ji = J(q(1), q(2), q(3), q(4), q(5));
        si_values = svd(Ji);
        val = min(si_values);
        sigmax = max(si_values);
        if val < 0
            val = 0;
        end
        kin_cost = 1 - (val/sigmax);
        a = (qu - ql)/2;
        b = (qu + ql)/2;
        bound_cost = (q - b)./a;
        values = kin_cost;
        for x = bound_cost'
            if x >= 0
                k = x;
                if x > 1
                    k = 1;
                end
                values = [values k];
            else
                k = x;
                if x < -1
                    k = -1;
                end
                values = [values -1*k];
            end
        end
        sing_val = -1 * max(values);
    end
    sing_func = @evaluate_singular_value;
end

function [nl_sing_con_func] = get_sing_nlcon(p, pi, val, r, p_old)
    function [c, ceq] = nlcon(q)
        pval = p(q(1), q(2), q(3), q(4), q(5));
        ceq = pval(pi) - val;
        dist_to_old = abs(norm(pval - p_old));
        c = [];
        c = dist_to_old - r;
    end
    nl_sing_con_func = @nlcon;
end