function [dA, ta] = generate_rand_dA_balanced()
    % return a random dA of a random taste, sampled uniformly or with the 
    % provided weight.
    
    ta = weighted_dice([1, 2, 3, 4, 5, 6], [0.2, 0.2, 0.2, 0.2, 0.2, 0.2]);
    % (sum of weights can be not normalized)
    dA = generate_rand_dA_by_taste(ta);
    
end
