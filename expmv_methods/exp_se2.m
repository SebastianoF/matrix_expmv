function A = exp_se2(dA)
    % function A = exp_se2_a(dA)
    % From a matrix of se2_a, or a vector with its parameters,
    % its return the closed form of its exponentiation.
    
    if size(dA, 1) == 3 && size(dA, 2) == 3
        theta = dA(2,1);
        dtx = dA(1,3);
        dty = dA(2,3);
    elseif size(dA, 1) == 1 && size(dA, 2) == 3
        theta = dA(1);
        dtx = dA(2);
        dty = dA(3);
    end

    if abs(theta) < 100 * eps
        tx = dtx;
        ty = dty;
    else
        f = 1/theta;
        tx = f * ( sin(theta)*dtx - (1 - cos(theta))* dty );
        ty = f * ( (1 - cos(theta))*dtx + sin(theta)* dty ); 
    end
    
    A = [cos(theta), -sin(theta), tx; ...
         sin(theta),  cos(theta), ty; ...
         0,           0,           1];
     
end