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

% CONSTRAINT 1:
% generate H in the Lie group so that its logm() has real entries 
% and h(d,d)=1.

Q = orth(rand(d));

v = abs(1+randn(1,d));  % diagonal vector
mv = max(v);
v(d) =(10*mv -1*(Q(d, 1:d-1).^2)*(v(1:d-1)'))/(Q(d,d)^2);
H = Q*diag(v)*Q';
H = H/(10*mv);

% or work with SPD:
% Q = rand(d);
% H = Q'*Q;

% centered H: according to the formula you provided me in the email
center = [10, 10]';

den = 1 - H(d, 1:d-1)*center;

A = (H(1:d-1, 1:d-1) + kron(center', H(d, 1:d-1)') ) / den;
B = H(d, 1:d-1) / den;
T = (-H(1:d-1, 1:d-1)*center - (B*center)*center + center) / den;

H_c = ones(d,d);
H_c(1:d-1, 1:d-1) = A;
H_c(d, 1:d-1) = B;
H_c(1:d-1, d) = T;

%H_c(d,d) = 1/den;

% Centered H should have the second CONSTRAINT: the real eigenvalues must 
% be positive.

h = logm(H);
h_c = logm(H_c);

% show matrices:

disp('h = ') 
disp(h)

disp('H = ') 
disp(H)

disp('h_c = ') 
disp(h_c)

disp('H_c = ') 
disp(H_c)

disp('norm(H-eye(d)), norm(H_c - eye(d)): ')
disp([norm(H-eye(d)), norm(H_c - eye(d)) ])

format short

    function dx = hom_vel(t,x)
        dx = zeros(size(x));
        dx(1,:) = h(1,1)*x(1,:) + h(1,2)*x(2,:) + h(1,3) - ...
                  x(1,:).*( h(3,1)*x(1,:) + h(3,2)*x(2,:) + h(3,3) );
        dx(2,:) = h(2,1)*x(1,:) + h(2,2)*x(2,:) + h(2,3) - ...
                  x(2,:).*( h(3,1)*x(1,:) + h(3,2)*x(2,:) + h(3,3) );
    end

end