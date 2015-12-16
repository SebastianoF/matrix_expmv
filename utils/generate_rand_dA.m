function [dA, c] = generate_rand_dA()
    
    % A = generate_rand_A()
    % generate a matrix of the form
    %
    % dA = [a, b, alpha;
    %       c, d, beta;
    %       0, 0, 0]
    %
    % with relative class c
    % according to the subdivision in classes provided in main.m 
    
    a = - 3;
    b = 3;
    S = (b - a)*rand(2, 2) + a;
    t = 5*rand(2,1);
    
    dA = [S,t];
    dA = [dA; 0, 0, 0];
    
    % Select the class of the generated matrix
    
    e = eig(S); l1 = e(1); l2 = e(2);
    
    % eig both reals
    if imag(l1) < (10^4)*eps && imag(l2) < (10^4)*eps
        l1 = real(l1); l2 = real(l2);
        
        % both real and same sings
        if l1*l2 > 0 
            c = 1;
        
            % both real and opposite sings
        else
            c = 2;
        end
        
    % eig complex conjugate
    else
        
        % real part is negative or positive
        if abs(real(l1)) > (10^4)*eps
            c = 3; 
            
        % real part is zero: 
        % NOTE this case will almost never be taken into account for 
        % random generated matrices!
        % real part of the eigenvalue are almost never equals to zero!
        % To generate a matrix of this type is necessarily to use the
        % function generate_se2 instead! See generate_rand_A_cl
        else
            c = 4;   
        end
    end
      
end