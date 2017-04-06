function hom_svf_2()

clc
format long

%%% options %%%

maxpix = 20;  % max x and max y domain of the svf
% projective dimension:
d = 3;
scale_factor = 1/maxpix;

% intial codition
x0 = 20; y0 = 1;

% Generate H as rotation, plus projective
theta = (pi/8)*rand(1);
tx = 0.5 * rand(1);
ty = 0.5 * rand(1);
pr1 = 0.5 * rand(1);
pr2 = 0.5 * rand(1);

theta_1 = (pi/12) * rand(1);

alpha = cos(theta_1);
beta = sin(theta_1);
delta = cos(theta_1);
gamma = sin(theta_1);

H = [alpha*cos(theta), -beta*sin(theta), tx; gamma*sin(theta), delta*cos(theta), ty; pr1, pr2, 1];

% Other constraint:

h = rand(d);

H = expm(h);


% show matrices:

disp('h = ') 
disp(h)

disp('H = ') 
disp(H)

disp('norm(H-eye(d)): ')
disp([norm(H-eye(d)) ])

format short

    function dx = hom_vel(t,x)
        dx = zeros(size(x));
        dx(1,:) = h(1,1)*x(1,:) + h(1,2)*x(2,:) + h(1,3) - ...
                  x(1,:).*( h(3,1)*x(1,:) + h(3,2)*x(2,:) + h(3,3) );
        dx(2,:) = h(2,1)*x(1,:) + h(2,2)*x(2,:) + h(2,3) - ...
                  x(2,:).*( h(3,1)*x(1,:) + h(3,2)*x(2,:) + h(3,3) );
    end

end