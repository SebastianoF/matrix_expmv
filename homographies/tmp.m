clc

a = rand(3,3);
ata = 100*rand(1)*a'*a;

h = randn(4);
h(4, 1:4) = abs(h(4, 1:4));

H = expm(h);


% show matrices:

disp('h = ') 
disp(h)

disp('H = ') 
disp(H)