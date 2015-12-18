function dA = generate_rand_dA_by_taste(taste_input)
    % Given the taste integer from 1 to 6 it generates the matrix of the 
    % shape of the appropriate class with a Montecarlo algoritm based on
    % the function generate_rand_A
    
    if taste_input == 6
        dA = generate_se2_dA();
    else

        flag = 0;
        while flag == 0
            [dA, ta] = generate_rand_dA();
            if ta == taste_input
                flag = 1;
            end
        end
        
    end
end