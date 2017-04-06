function homsvf_test()

A1 = [0 0 1; 
      0 0 0; 
      0 0 0];

A2 = [0 0 0; 
      0 0 1; 
      0 0 0];

A3 = [0 1 0; 
      0 0 0; 
      0 0 0];

A4 = [0 0 0; 
      1 0 0; 
      0 0 0];

A5 = [1 0 0; 
      0 -1 0; 
      0 0 0];

A6 = [0 0 0; 
      0 -1 0; 
      0 0 1];
  
A7 = [0 0 0; 
      0 0 0; 
      1 0 0];

A8 = [0 0 0; 
      0 0 0; 
      0 1 0];

a = [10; 10; 0.2; 0.2; 0.2; 0.2; 0.05; 0.05].*randn(8,1)*10;
% a(1) = 10;
% a(2) = -5;
% a(3) = 0.2;
% a(4) = -0.2;
% a(5) = 0.3;
% a(6) = 2;
% a(7) = -0.2;
% a(8) = 0.01;

A = a(1)*A1 + a(2)*A2 + a(3)*A3 + a(4)*A4 + a(5)*A5 + a(6)*A6 + a(7)*A7 + a(8)*A8

options = odeset('RelTol',1e-4,'AbsTol',[1e-4 1e-4]);

x0 = 5*randn(2,1);
%x0 = [2; 10]

% computed solution:
[T,X] = ode45(@hom_vel, [0 1], x0, options);
x1 = X(end,:)'

% analytic solution
eaprod = expm(A)*[x0; 1];

x1_analytic = eaprod(1:2)./eaprod(3)

maxpix = 10;
[xx,yy] = meshgrid(0:maxpix,0:maxpix);
svf = hom_vel(0,[xx(:)'; yy(:)']);
svfx = reshape(svf(1,:),maxpix+1,maxpix+1);
svfy = reshape(svf(2,:),maxpix+1,maxpix+1);

quiver(xx,yy,svfx,svfy,0)

function dx = hom_vel(t,x)
dx = zeros(size(x));
dx(1,:) = a(1) + a(3)*x(2,:) + a(5)*x(1,:) ...
    - x(1,:).*( a(6) + a(7)*x(1,:) + a(8)*x(2,:) );
dx(2,:) = a(2) + a(4)*x(1,:) - a(5)*x(2,:) - a(6)*x(2,:) ...
    - x(2,:).*( a(6) + a(7)*x(1,:) + a(8)*x(2,:) );
end

end