function face = weighted_dice(faces, weight)
    % face = dice(faces, weight) 
    % input: row vector. Faces and weight of the same dimension!
    % output: face is the discrete face after tossing a virtual dice with given
    % faces and weight.
    % thanks to http://codereview.stackexchange.com/questions/110915/weighted-dice-with-n-faces/111030#111030
    
    if size(faces, 2) ~= size(weight, 2) || size(weight, 2)<= 1
        error('Input of function dice not well defined. See help')
    end

    num_weight = size(weight, 2);

    cumulative_weight = zeros(1, num_weight + 1);
    for j=1:num_weight
       cumulative_weight(j+1) = sum(weight(1:j));
    end
    rand_w = unifrnd(0, cumulative_weight(end));
    for j=1:num_weight
       if rand_w >= cumulative_weight(j) && rand_w < cumulative_weight(j+1)
           break
       end
    end

    face = faces(j);