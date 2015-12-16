
clear 
%disp('One answer: ')
%disp(weighted_dice(['a','b','c','d','e','f'], [0.1, 0.1, 0.5, 0.1, 0.1, 0.1]))

beans = zeros(1,6);
for i=1:1000
    val = weighted_dice(['a','b','c','d','e','f'], [0.1, 0.1, 0.5, 0.1, 0.1, 0.1]);
    if val == 'a'
        beans(1,1) = beans(1,1) + 1;
    elseif val == 'b'
        beans(1,2) = beans(1,2) + 1;
    elseif val == 'c'
        beans(1,3) = beans(1,3) + 1;
    elseif val == 'd'
        beans(1,4) = beans(1,4) + 1;
    elseif val == 'e'
        beans(1,5) = beans(1,5) + 1;
    elseif val == 'f'
        beans(1,6) = beans(1,6) + 1;
    end
end

%disp(beans)

figure(4);
bar(beans)
set(gca,'XTickLabel',{'0.1', '0.1', '0.5', '0.1', '0.1', '0.1'});
xlabel('weight of each face');
ylabel('sampling frequency');

if beans(3) < 400
    disp('test_dice NOT passed, or very bad luck (?)');
else
    disp('test_dice 0 passed');
end




