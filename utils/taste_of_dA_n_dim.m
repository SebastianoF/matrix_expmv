function ta = taste_of_dA_n_dim(dA)
    
    % ta = taste_of_dA_n_dim(dA)
    % Returns the taste of dA with the multidimensional nomenclature.
    %
    % Taste 'neg'
    % If all the eigenvalues, both real or complex, have negative or zero 
    % real part. (zero real part is unlikely to happen)
    % 
    % Taste 'pos'
    % If all the eigenvalues, both real or complex, have positive real part. 
    % 
    % Taste 'mix'
    % If all the eigenvalues, both real or complex, have bot positive and 
    % negative real part.
    % 

    assert(size(dA, 1) == size(dA, 2))
    
    e = eig(dA);
    all_neg = true;
    all_pos = true;
    
    for j=1:size(e,1)
       if real(e(j)) > 0
           all_neg = false;
       elseif real(e(j)) < 0
           all_pos = false;
       end
    end
    
    if all_neg == false && all_pos == false
        ta = 'mix';
    elseif all_neg == true && all_pos == false
        ta = 'neg';
    elseif all_neg == false && all_pos == true
        ta = 'pos';
    else
        % real part are all zero
        ta = 'neg';
    end

end