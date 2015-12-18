function dA_large = generate_rand_dA_by_taste_n_dim(n, taste_input)
    % Given a dimension and a taste it generates a nxn random matrix
    % of appropriate taste. See taste nomenclature for big matrices on the
    % README.md file.
    % For the moment only real eigenvalues random matrices are created.

    if strcmp('pos', taste_input)
        D = 10 * diag(abs(randn(n,1))); 
        P = orth(randn(n)); 
        dA_large = P*D*P';
    elseif strcmp('neg', taste_input)
        D = 10 * diag(-abs(randn(n,1))); 
        P = orth(randn(n)); 
        dA_large = P*D*P';
    elseif strcmp('mix', taste_input)
        a = datasample([1,-1],n) .* randn(1,n);
        D = 10 * diag(a); 
        P = orth(randn(n)); 
        dA_large = P*D*P';
    end
    
end