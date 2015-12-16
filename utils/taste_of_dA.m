function ta = taste_of_dA(dA)
    % Given a matrix that defines 2d linear stationary ode system in 
    % homogeneous coordinates, returns its taste, according to the
    % nomenclature (see README.dm of the project).
    
    S = dA(1:2, 1:2);
    
    e = eig(S); l1 = e(1); l2 = e(2);
    
    % eig both reals
    if imag(l1) < (10^4)*eps && imag(l2) < (10^4)*eps
        l1 = real(l1); l2 = real(l2);
        
        % both real and positive
        if l1 > 0 && l2 > 0 
            ta = 1;
        
        % both real and negative
        elseif l1 < 0 && l2 < 0 
            ta = 2;
        
            % both real opposite signs
        else
            ta = 3;
            
        end
        
    % eig complex conjugate
    else
        
        % real part is positive
        if real(l1) > (10^4)*eps
            ta = 4; 
            
        % real part is negative
        elseif real(l1) < -(10^4)*eps
            ta = 5; 
            
        % real part is zero: 
        % NOTE this case will almost never be taken into account for 
        % random generated matrices!
        % real part of the eigenvalue are almost never equals to zero!
        % To generate a matrix of this type is necessarily to use the
        % function generate_se2 instead! See generate_rand_A_by_taste.
        else
            ta = 6;   
        end
    end
end