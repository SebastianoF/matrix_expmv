function homsvf_test()

format long

%%% options %%%

maxpix = 20;  % max x and max y domain of the svf
scale_factor = 1/maxpix;
% intial codition
x0 = 20; y0 = 1;

opt = 3;

%%% example 1:
if opt == 1 
    % good option obtained with an trick!
    h = (1 / maxpix) * 0.5 * randn(3,3);
    n = max(max(h));
    h(3,:) = [n n n];
    
elseif opt == 2
    % bad option
    h = scale_factor * [8.943142  2.182549  4.82487;
                        -9.317464 -1.386941 -1.773795;
                        -4.13707 -3.135003 -2.19091];
    
elseif opt == 3
    % graphically pleasant option, but it creates problems on the border
    % when compared with the scaling and squaring.
    % When the norm of the arrays peakes around the border we have problem.
    % the scaling and squaring moves much slower than the analytic
    % solution.
    
    h = scale_factor *10* [8.943142  -2.182549  4.82487;
                        -9.317464 1.386941 -1.773795;
                        4.13707 -3.135003 2.19091];
                    
elseif opt == 4
    % not good option (belonging to psl(3) does not guarantee better restults)
    h =  [0.5      0     0.5; 
          -0.5     0.5   -0.50; 
          -.5        0     -1.0];
      
elseif opt == 5
    % good option on the border
    x0 = 2; y0 = 19;
    h = [ 0.08943142    0.02182549    0.00482487;
      -0.09317464   -0.01386941   -0.01773795;
      -0.00413707   -0.03135003   -0.00219091];

else
    h = scale_factor * 0.5 * randn(3,3);
end
 
disp('h = ') 
disp(h)

H = expm(h);

disp('H = ') 
disp(H)

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

    function dx = hom_vel(t,x)
        dx = zeros(size(x));
        dx(1,:) = h(1,1)*x(1,:) + h(1,2)*x(2,:) + h(1,3) - ...
                  x(1,:).*( h(3,1)*x(1,:) + h(3,2)*x(2,:) + h(3,3) );
        dx(2,:) = h(2,1)*x(1,:) + h(2,2)*x(2,:) + h(2,3) - ...
                  x(2,:).*( h(3,1)*x(1,:) + h(3,2)*x(2,:) + h(3,3) );
    end

end