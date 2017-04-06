


function hom_svf_from_group()

clc
format long

%%% options %%%

maxpix = 20;  % max x and max y domain of the svf
% projective dimension:
d = 3;
scale_factor = 1/maxpix;

% intial codition
x0 = 20; y0 = 1;

% generate h in the Lie group (its logm() must have real entries 
% and h(2,2)=1). This restriction may cause problems...

Q = orth(rand(d));

v = abs(1+randn(1,d));  % diagonal vector
mv = max(v);
v(d) =(10*mv -1*(Q(d, 1:d-1).^2)*(v(1:d-1)'))/(Q(d,d)^2);
H = Q*diag(v)*Q';
H = H/(10*mv);





% TODO: from here. Visualize matrices 
% then do a second version where the SVF is centered.

% centered H:
center = [10, 10]';

den = 1 - H(d, 1:d-1)*center;

A = (H(1:d-1, 1:d-1) + kron(center', H(d, 1:d-1)') ) / den;
B = H(d, 1:d-1) / den;
T = (-H(1:d-1, 1:d-1)*center - (B*center)*center + center) / den;

H_c = ones(d,d);
H_c(1:d-1, 1:d-1) = A;
H_c(d, 1:d-1) = B;
H_c(1:d-1, d) = T;

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

% intial point
p0 = [x0, y0]';  % 5*randn(2,1);

msg = sprintf('svf at the point x,y = %d, %d', x0, y0);
disp(msg)
disp(hom_vel(1, [x0,y0,1]'))

msg = sprintf('displacement at the point x,y = %d, %d', x0, y0);
disp(msg)

% computed solution:
options = odeset('RelTol',1e-4,'AbsTol',[1e-4 1e-4]);
[T,X] = ode45(@hom_vel, [0 1], p0, options);

svf_at_p0 = hom_vel(1,[p0; 1]);

p1 = X(end,:)';

disp('p1 computed with ode45:')
disp(p1)

% analytic solution:
S = expm(h)*[p0; 1];

% disp('Analytic solution in projective coordinates S = expm(h)*[x0; 1]')
% disp(S)

p1_analytic = S(1:2)./S(3);

disp('analytic solution in affine coordinates S[0:2]/S[2] = :')
disp(p1_analytic)

% generate integral curve
times = 0:0.05:1;
num_time_points = size(times, 2);

integral_curve = zeros(3, num_time_points);

for t = 1:num_time_points
    integral_curve(:, t) = expm(times(t) * h) * [x0; y0; 1];
end

% generate SFV
[xx,yy] = meshgrid(0:maxpix,0:maxpix);
svf = hom_vel(0,[xx(:)'; yy(:)']);
svfx = reshape(scale_factor * svf(1,:),maxpix+1,maxpix+1);
svfy = reshape(scale_factor * svf(2,:),maxpix+1,maxpix+1);

% plot
figure(11)
quiver(xx, yy, svfx, svfy, 0, 'r')
hold on
% quiver(x0, y0, ...
%        scale_factor * p1_analytic(1), scale_factor * p1_analytic(2), 'm', 'linewidth', 3);
% plot(integral_curve(1, :), integral_curve(2, :), 'b', 'linewidth', 4)
hold off

format short

    function dx = hom_vel(t,x)
        dx = zeros(size(x));
        dx(1,:) = h(1,1)*x(1,:) + h(1,2)*x(2,:) + h(1,3) - ...
                  x(1,:).*( h(3,1)*x(1,:) + h(3,2)*x(2,:) + h(3,3) );
        dx(2,:) = h(2,1)*x(1,:) + h(2,2)*x(2,:) + h(2,3) - ...
                  x(2,:).*( h(3,1)*x(1,:) + h(3,2)*x(2,:) + h(3,3) );
    end

end