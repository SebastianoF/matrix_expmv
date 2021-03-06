function dA = generate_se2_dA(param)
    % A = generate_se2_a(theta, tx, ty)
    % generate an element of the Lie algebra se2
    % The corrseponding element in the Lie group can be generated using
    % a closed form.
    % If no input arguments are provided then the matrix is random,
    if nargin == 0
        a = - 5; 
        b = 5;
        theta = (b-a)*rand() + a;
        tx = 10*rand()-5;
        ty = 10*rand()-5;
    elseif nargin == 1 && size(param, 2) == 3 
        theta = param(1);
        tx = param(2);
        ty = param(3);
    else
        error('Input values are not correct')
    end
    dA = [0, -theta , tx; ...
          theta, 0,    ty; ...
          0,     0,    0];
   
end