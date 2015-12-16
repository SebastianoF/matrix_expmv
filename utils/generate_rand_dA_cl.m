function dA = generate_rand_dA_cl(c_input)
    % Given the class c it generates the matrix of the shape
    % of the appropriate class with a Montecarlo algoritm based on the
    % function generate_rand_A
    
    if c_input == 4
        dA = generate_se2_dA();
    else

        flag = 0;
        while flag == 0
            [dA, c] = generate_rand_dA();
            if c == c_input
                flag = 1;
            end
        end
        
    end
end