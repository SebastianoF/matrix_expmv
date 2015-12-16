function [dA, ty] = generate_rand_dA()
    
    % A = generate_rand_A()
    % generate a matrix of the form
    %
    % dA = [a, b, alpha;
    %       c, d, beta;
    %       0, 0, 0]
    %
    % with relative type ty
    % according to the subdivision in classes provided in main.m 
    
    a = - 5;
    b = 5;
    S = (b - a)*rand(2, 2) + a;
    t = 10*rand(2,1)-5;
    
    dA = [S,t];
    dA = [dA; 0, 0, 0];
    
    % Select the class of the generated matrix
    
    e = eig(S); l1 = e(1); l2 = e(2);
    
    % eig both reals
    if imag(l1) < (10^4)*eps && imag(l2) < (10^4)*eps
        l1 = real(l1); l2 = real(l2);
        
        % both real and positive
        if l1 > 0 && l2 > 0 
            ty = 1;
        
            % both real and negative
        elseif l1 < 0 && l2 < 0 
            
        else
            ty = 3;
        end
        
    % eig complex conjugate
    else
        
        % real part is negative or positive
        if abs(real(l1)) > (10^4)*eps
            ty = 3; 
            
        % real part is zero: 
        % NOTE this case will almost never be taken into account for 
        % random generated matrices!
        % real part of the eigenvalue are almost never equals to zero!
        % To generate a matrix of this type is necessarily to use the
        % function generate_se2 instead! See generate_rand_A_cl
        else
            ty = 4;   
        end
    end
      
end