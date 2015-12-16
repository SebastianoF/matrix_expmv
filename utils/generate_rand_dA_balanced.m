function [dA, c] = generate_rand_dA_balanced()
    % return a random dA of a random class, sampled uniformly or with the 
    % provided weight.
    
    c = weighted_dice([1, 2, 3, 4], [0.25, 0.25, 0.25, 0.25]);
    dA = generate_rand_dA_cl(c);
    
end
