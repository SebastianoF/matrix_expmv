function [dA, ta] = generate_rand_dA()
    
    % A = generate_rand_A()
    % generate a matrix of the form
    %
    % dA = [a, b, alpha;
    %       c, d, beta;
    %       0, 0, 0]
    %
    % with output taste ty
    % according to the subdivision in tastes provided in README.m 
    
    a = - 5;
    b = 5;
    S = (b - a)*rand(2, 2) + a;
    t = 10*rand(2,1)-5;
    
    dA = horzcat(S,t);
    dA = vertcat(dA, [0, 0, 0]);
    
    % Return the taste of the generated matrix
    ta = taste_of_dA(dA);
      
end