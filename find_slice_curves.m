function [xy,Theta_min, Theta_max] = find_slice_curves(num_of_scurves, major_search, minor_search, T, ll) % Arguments?
    % we will do the calculations via symbolic objects & inline functions
    x = matlabFunction(T(1,4)); % depends on (th1,th2,th3,th4,th5)
    y = matlabFunction(T(2,4)); % depends on (th1,th2,th3,th4,th5)
    z = matlabFunction(T(3,4)); % depends on (th2,th3,th4,th5)
    
    vars_funcs = {x, y, z};
    coord = ['x' 'y' 'z'];
    % getting the objective functions for minimizing & maximizing
    calculate_x_val = get_objective_function(1, vars_funcs, major_search);
    calculate_x_val_maximizing = get_objective_function(-1, vars_funcs, major_search);
    calculate_y_val = get_objective_function(1, vars_funcs, minor_search);
    calculate_y_val_maximizing = get_objective_function(-1, vars_funcs, minor_search);

    % here we carry out the actual optimization.
    % Although fmincon() is designed for minimizing, it can also be used for maximizing
    % as well by modifying the objective function (negating it). This is
    % done through get_objective_func
    opts = optimset('Display','off','MaxIter', 3000, 'MaxFunEvals', 3000);
    
    Th_min = fmincon(calculate_x_val,ll(1:5,1)',[],[],[],[],ll(1:5,1)',ll(1:5,2)', [], opts);   
    Th_max = fmincon(calculate_x_val_maximizing,ll(1:5,1)',[],[],[],[],ll(1:5,1)',ll(1:5,2)', [], opts);
    
    % calculate the minimum and maximum values
    x_min = calculate_x_val(Th_min);
    x_max = calculate_x_val(Th_max);
    
    fprintf(['The minimum ' coord(major_search) ' config is at (%0.2f, %0.2f, %0.2f, %0.2f, %0.2f) where' coord(major_search) ' = %0.2f\n'], [Th_min x_min]);
    fprintf(['The maximum ' coord(major_search) ' config is at (%0.2f, %0.2f, %0.2f, %0.2f, %0.2f) where' coord(major_search) ' = %0.2f\n'], [Th_max x_max]);
    
    Theta_min = [];
    Theta_max = [];
    xy = [];
    counter = 1;
    for xv = linspace(x_min, x_max, num_of_scurves)
        x_nlcon = get_nlcon_function(xv, vars_funcs, major_search);
        Th_min = fmincon(calculate_y_val,ll(1:5,1)',[],[],[],[],ll(1:5,1)',ll(1:5,2)', x_nlcon, opts);   
        Th_max = fmincon(calculate_y_val_maximizing,ll(1:5,1)',[],[],[],[],ll(1:5,1)',ll(1:5,2)', x_nlcon, opts);
        Theta_min = [Theta_min; Th_min];
        Theta_max = [Theta_max; Th_max];
        y_min = calculate_y_val(Th_min);
        y_max = calculate_y_val(Th_max);
        xy = [xy; xv y_min y_max];
        fprintf('%d : (min, max) %c value is (%0.2f, %0.2f) at %c = %0.2f\n',counter, coord(minor_search), y_min, y_max, coord(major_search), xv);
        counter = counter + 1;
    end
end

function [xyz_func] = get_objective_function(sign, vars_funcs, func_index)
    function [xyz_val] = calculate_xyz_val(theta_val)
        var_func = vars_funcs{func_index};
       
        t1 = theta_val(1);
        t2 = theta_val(2);
        t3 = theta_val(3);
        t4 = theta_val(4);
        t5 = theta_val(5);
        if func_index == 3
            xyz_val = sign * var_func(t2, t3, t4, t5);
        else
            xyz_val = sign * var_func(t1, t2, t3, t4, t5);
        end
      
    end

    xyz_func = @calculate_xyz_val;
end

function [nlcon] = get_nlcon_function(cval, vars_funcs, func_index)
    function [c, ceq] = nlcon_func(theta_val)
        c = [];
        var_func = vars_funcs{func_index};
        t1 = theta_val(1);
        t2 = theta_val(2);
        t3 = theta_val(3);
        t4 = theta_val(4);
        t5 = theta_val(5);
        if func_index == 3
            ceq = cval - var_func(t2, t3, t4, t5);
        else
            ceq = cval - var_func(t1, t2, t3, t4, t5);
        end
    end
    nlcon = @nlcon_func;
end