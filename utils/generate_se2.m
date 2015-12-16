function dA = generate_se2(param)
    % A = generate_se2_a(theta, tx, ty)
    % generate an element of the Lie algebra se2
    % The corrseponding element in the Lie group can be generated using
    % a closed form.
    % If no input arguments are provided then the matrix is random,
    if nargin == 0
        a = - pi/4; 
        b = pi/4;
        theta = (b-a)*rand() + a;
        tx = 5*rand();
        ty = 5*rand();
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