function [phikAv, errest, info, extreigs, inieigs] = phileja (h, A, varargin)
% [PHIKAV, ERREST, INFO, EXTREIGS] = PHILEJA (H, A, K, V, TOL, EXTREIGS)
%
% Compute phi_K(H * A) * V, where A is a matrix *preferably* sparse and with
% eigenvalues with negative real part. V is a vector or a matrix.
% If length (TOL) == 1, then TOL is the absolute tolerance. If 
% length (TOL) == 2, then TOL(1) is the absolute tolerance and TOL(2)
% the tolerance relative to the current approximation of phi_K(H * A)*V. If
% not given, then TOL = [2^(-24), 0].
%
% NEXT FEATURE NOT WORKING YET (July 14, 2011)
%
% H can be a row vector: in this case, the result is equivalent
% to [PHILEJA (H(1), A, K, V, ...),...,PHILEJA (H(end), A, K, V, ...)],
% but much faster. If H is a row vector and V a matrix with 
% size(V, 2) == length (H), the result is
% [PHILEJA (H(1), A, K, V(:,1), ...),...,PHILEJA (H(end), A, K, V(:,end), ...)].
%
% END
%
% If EXTREIGS is not given, the spectrum of A is estimated by
% Gershgorin's disks. Otherwise, it is estimated by EXTREIGS.SR
% (smaller real part of the eigenvalues), EXTREIGS.LR (largest real
% part) and EXTREIGS.LI2 (squared largest imaginary part).
% On output, sum (INFO) is the total number of iterations 
% (== matrix-vector products) and sum(ERREST) the estimated error (in
% infinity norm). On output, EXTREIGS is like above.
%
% [PHIKAV, ERREST, INFO, EXTREIGS, INIEIGS] = PHILEJA (H, AF, ISSYM, K, V, TOL)
%
% Compute phi_K(H * A)(V), where AF is a function handle for the linear
% operator A(V). ISSYM is logical and specifies if AF is symmetric or not.
% On output, EXTREIGS is like above and INIEIGS is a structure with  
% vectors that can be used in a next call to PHILEJA to better estimate 
% extreme eigenvalues. If AF is symmetric, INIEIGS is a structure with the 
% unique field BE (vector of size SIZE (V)). If AF is not symmetric, it is a 
% structure with the fields SR, LR and LI (all of them vectors of size
% SIZE (V)). If EXTREIGS or INIEIGS are already known structures, it is 
% possible to call 
%
% PHILEJA (H, AF, ISSYM, K, V, TOL, EXTREIGS)
%
% (in this case, on output INIEIGS is empty) or
%
% PHILEJA (H, AF, ISSYM, K, V, TOL, INIEIGS)
%
% [PHIKAV, ERREST, INFO, EXTREIGS, INIEIGS] = PHILEJA (H, AFSYM, AFSSYM, K, V, TOL)
%
% Compute phi_K(H* (AFSYM + AFSSYM))(V), where AFSYM is a function handle for
% the symmetric part of the linear operator A(V) and AFSSYM is a
% function handle for the skew-symmetric part of A(V). On output,
% EXTREIGS is like above and INIEIGS is a structure with vectors that
% can be used in a next call to PHILEJA to better estimate extreme 
% eigenvalues. It has the fields BE and LM. If EXTREIGS or INIEIGS are
% already known structures, it is possible to call 
%
% PHILEJA (AFSYM, AFSSYM, K, V, TOL, EXTREIGS)
%
% (in this case, on output INIEIGS is empty) or
%
% PHILEJA (AFSYM, AFSSYM, K, V, TOL, INIEIGS)

% spectrum estimate
if (isfloat (A))
  v = varargin{2};
  if (length (v) ~= size (A, 2))
    size(v),size(A)
    error ('Inconsistent matrix and vector sizes')
  end
  k = varargin{1};
  if (nargin >= 5)
    tol = varargin{3};
    if (length (tol) == 1)
      tol(2) = 0;
    end
  else
    tol = [2^(-24),0];
  end
  if (nargin == 6)
    extreigs = varargin{4};
  else
% A is a matrix, use Gershgorin's disks
    extreigs = gersh(A);
  end
  inieigs = [];
  A = @(u) (A * u) * diag (h);
else
  if (nargin < 5 || ...
      (~isa (varargin{1}, 'logical') && ~isa (varargin{1}, 'function_handle')))
    error ('call phileja (AFsym, AFskewsym, k, v) or phileja (AFsym, [], k, v)')
  elseif (isa (varargin{1}, 'function_handle'))
% symmetric and skewsymmetric part given
    k = varargin{2};
    v = varargin{3};
    n = length(v);
    AFskewsym = varargin{1};
    if (nargin == 7)
      if (isfield (varargin{5}, 'LI2'))
% eigenvalues
	extreigs = varargin{5};
	inieigs = [];
      else
% initial eigenvectors	
	inieigs = varargin{5};
	[extreigs, inieigs] = eigssymskewsym (n, A, AFskewsym, inieigs);
      end	
    else  
      [extreigs, inieigs] = eigssymskewsym (n, A, AFskewsym);
    end  
    A = @(u) A (u) * diag (h) + AFskewsym (u) * diag (h);
  else
% only a single linear operator
    k = varargin{2};
    v = varargin{3};
    n = length(v);
    if (nargin == 7)
      if (isfield (varargin{5}, 'LI2'))
% eigenvalues
	extreigs = varargin{5};
	inieigs = [];
      else
% initial eigenvectors	      
	inieigs = varargin{5};
	[extreigs, inieigs] = eigsfull (n, A, varargin{1}, inieigs);
      end	
    else      
      [extreigs, inieigs] = eigsfull (n, A, varargin{1});
    end  
    A = @(u) A(u) * diag(h);
  end
  if (nargin >= 6 && ~isempty(varargin{4}))
    tol = varargin{4};
    if (length (tol) == 1)
      tol(2) = 0;
    end
  else  
    tol = [2^(-24),0];
  end
end
SR = extreigs.SR * max(h);
LR = extreigs.LR * max(h);
LI2 = extreigs.LI2 * max(h) ^ 2;
d = (SR + LR) / 2;
C2 = ((LR - SR) / 2) ^ 2;
c2 = (nthroot (C2, 3) + nthroot (LI2, 3)) * (C2 ^ (2/3) - LI2 ^ (2/3));
gamma = (nthroot (C2, 3) * sqrt (nthroot (C2, 3) + nthroot (LI2, 3)) + ...
    sqrt (nthroot (C2 * LI2 ^ 2, 3) + LI2)) / 2;
% inscribed ellipse
%a = -SR/2;
%b = sqrt(LI2);
%gamma = (a+b)/2
%c2 = a^2-b^2
% in order to treat trivial cases gammafoc = 0 (diagonal matrices)
gammafoc = max (sqrt (abs (c2)) / 2, 1);
%gammafoc = max(-d/2,1) % to interpolate exactly h=0
if (gamma > 10 * gammafoc)
% circle case
%  gammafoc = min (gamma, sqrt (abs (c2)) / 2);
  gammafoc = min (gamma, abs(d) / 2);
end
% three magic numbers...
maxm = 150;
minm = 30;
m = max (2.6 * ceil (gamma), minm); 
% end
nsteps = ceil (m / maxm);
tol = tol / nsteps;
hinner = 1 / nsteps;
if (length (h) > 1)
  if (size (v, 2) == 1)
    v = repmat (v, 1, length (h));
  elseif (size (v, 2) ~= length (h))
    error('size(V, 2) must be equal to length (H)')
  end
end
m = maxm;
%gamma
%m = 40;
%nsteps = 1;
%hinner = 1;
if (c2 >= 0)
% real Leja points
  xi = leja(m) + d / gammafoc;
  newt = @newton;
else
%  disp('complex')
% symmetriezed complex Leja points
  xi = 1i * lejas(m) + d / gammafoc;
%  plot(gammafoc*real(xi),gammafoc*imag(xi),'o')
  newt = @newtons;
%  xi = 1i * leja(m) + d / gammafoc;
%  newt = @newton;
end
if (LI2 > 0)
% Fast complex Leja points
%  disp('fast')
%  newt = @newton;
%  gammafoc = gamma;
%  xi = flp (d,c2,gamma,m) + d / gamma;
%figure
%plot(gamma*real(xi),gamma*imag(xi),'*')
%stop
end
%dd = divdif(@phiktaylor,k,xi,hinner,gammafoc);
%maxm = find(dd < 0 | abs (dd) < eps,1);
%m = max (5 * ceil (gamma), minm); % magic number
%nsteps = ceil (m / maxm);
hinner = 1 / nsteps;
m = floor (m / nsteps);
%gammafoc = min (- d / 2, gammafoc);
%c = sqrt(c2);
%LI = sqrt(LI2);
%a = gamma + c2/(4 * gamma);
%b = gamma - c2/(4 * gamma);
%plot(SR,LI,'b*',LR,LI,'b*',SR,-LI,'b*',LR,-LI,'b*')
%hold on
%plot(d+a,0,'ro',d-a,0,'ro',d,b,'ro',d,-b,'ro')
%plot(d-c,0,'g+',d+c,0,'g+')
%w = exp(1i*linspace(0,2*pi,100));
%z = gamma*(w+d/gamma+c^2/(4*gamma^2)./w);
%plot(real(z),imag(z))
%z = gammafoc*(w+d/gammafoc+c^2/(4*gammafoc^2)./w);
%plot(real(z),imag(z),'*')
% usare time marching
%dd = divdif(@phiktaylor,0,xi,hinner,gammafoc)
%[phikAv,errest,info] = newt(A,v,xi,dd,gammafoc,tol(2),tol(1)); 
%info
%dd = divdif(@phiktaylor,1,xi,hinner,gammafoc)
%[phikAv,errest,info] = newt(A,v,xi,dd,gammafoc,tol(2),tol(1)); 
%info
%stop
if (k == 0)
  dd = divdif(@phiktaylordd,k,xi,hinner,gammafoc);
  phikAv = v;
  errest = [];
  info = [];
  for j = 1:nsteps
    [phikAv,perrest,pinfo] = newt(A,phikAv,xi,dd,gammafoc,tol(2),tol(1)); 
    errest = [errest,perrest];
    info = [info, pinfo];
  end   
elseif (k > 0)
  j = 1;
  dd = divdif(@phiktaylordd,k,xi,hinner,gammafoc);
  errest = [];
  info = [];
  [y{k},perrest,pinfo] = newt(A,v * diag (hinner .^ k),...
			       xi,dd(:,k+1),gammafoc,tol(2),tol(1)); 
  errest = [errest,perrest];
  info = [info, pinfo];
  phikAv = y{k};
  if (nsteps > 1) 
    j = 2;
    t = (j - 1) * hinner;
    for i = k - 1:-1:1
%      dd = divdif(@phiktaylor,i,xi,hinner,gammafoc);
      [y{i},perrest,pinfo] = newt(A,v / factorial(k-i) * diag (hinner .^ i),...
				   xi,dd(:,k+1),gammafoc,tol(2),tol(1)); 
      errest = [errest,perrest];
      info = [info, pinfo];
    end	
%    dd = divdif(@phiktaylor,0,xi,hinner,gammafoc);
    [phikAv,perrest,pinfo] = newt(A,phikAv,...
				 xi,dd(:,1),gammafoc,tol(2),tol(1)); 
    errest = [errest,perrest];
    info = [info, pinfo];
    for i = 1:k
      y{i} = y{i} * t ^ (k - i);
    end
%    phikAv = phikAv + sum(y * diag(t.^[k-1:-1:0]),2);
    phikAv = phikAv + sum (cat (3, y{:}), 3);
    for j = 3:nsteps
      [phikAv,perrest,pinfo] = newt(A,phikAv,...
				   xi,dd(:,1),gammafoc,tol(2),tol(1)); 
      errest = [errest,perrest];
      info = [info, pinfo];
      for i = 1:k
	y{i} = y{i} * ((j - 1) / (j - 2)) ^ (k - i);
      end
      phikAv = phikAv + sum (cat (3, y{:}), 3);
    end
  end  
else
  error(sprintf('phi_k for k = %d not defined',k))
end

function extreigs = gersh (A)
% The matrix is split into the symmetric and the skew-symmetric part
% Copyright Marco Caliari <marco.caliari@univr.it> 2010
% LR, SR (real part), LI (imaginary part).
AS = (A + A') / 2;
radius = sum (abs (AS), 2);
extreigs.SR = full (min (diag (AS) - (radius - diag (abs (AS)))));
extreigs.LR = min (full (max (diag (AS) + (radius - diag (abs (AS))))), 0);
AS = A - AS;
LI = full (max (diag (AS) + (sum (abs (AS), 2) - diag (abs (AS)))));
extreigs.LI2 = LI ^ 2;

function [extreigs, inieigs] = eigsfull (n, A, symmetry, inieigs)
% use eigs  
opts.maxit = 1;
opts.tol = 1;
opts.isreal = true;
opts.issym = symmetry;
opts.disp = 0;
if (symmetry)
  if (nargin == 4)
    opts.v0 = inieigs.BE;
  else
    opts.v0 = linspace(1,n,n)';
  end  
  [V, D, flag] = eigs (A, n, min (10, n - 2), 'BE', opts);
  eigenvalues = diag (D);
  [extreigs.SR, ind] = min (eigenvalues);
  inieigs.BE = V(:,ind);
  extreigs.LR = min (max (eigenvalues), 0);
  extreigs.LI2 = 0;
else
  if (nargin == 4)
    opts.v0 = inieigs.SR;
  else
    opts.v0 = linspace(1,n,n)';
  end  
  [V, D, flag] = eigs (A, n, min (10,n-2), 'SR', opts);
  eigenvalues = diag (D);
  [extreigs.SR, ind] = min (real (eigenvalues));
  inieigs.SR = V(:,ind);
  if (nargin == 4)
    opts.v0 = inieigs.LR;
  else
    opts.v0 = linspace(1,n,n)';
  end  
  [V, D, flag] = eigs (A, n, min (10,n-2), 'LR', opts);
  eigenvalues = diag (D);
  [LR, ind] = max (real (eigenvalues));
  extreigs.LR = min (LR, 0);
  inieigs.LR = V(:,ind);
  if (nargin == 4)
    opts.v0 = inieigs.LI;
  else
    opts.v0 = linspace(1,n,n)';
  end  
  [V, D, flag] = eigs (A, n, min (10, n - 2), 'LI', opts);
  eigenvalues = diag (D);
  [LI2, ind] = max (imag (eigenvalues));
  extreigs.LI2 = LI2 ^ 2;
% eigs requires a real initial vector for a real problem
  inieigs.LI =  real(V(:,ind));
end

function [extreigs, inieigs] = eigssymskewsym (n, A, AFskewsym, inieigs)
% use eigs on the symmetric part and the skewsymmetric part
% extreigs0 on input contains initial eigenvectors
opts.maxit = 10;
opts.tol = 1;
opts.isreal = true;
opts.issym = true;
opts.disp = 0;
if (nargin == 4)
  opts.v0 = inieigs.BE;
else
  opts.v0 = linspace(1,n,n)';
end
AF = @(u) A (u);
[V, D, flag] = eigs (AF, n, min (10, n - 2), 'BE', opts);
eigenvalues = diag (D);
[extreigs.SR, ind] = min (eigenvalues);
inieigs.BE = V(:, ind);
extreigs.LR = min (max (eigenvalues), 0);
AF = @(u) AFskewsym (AFskewsym (u));
if (nargin == 4)
  opts.v0 = inieigs.LM;
else
  opts.v0 = linspace(1,n,n)';
end
[V, D, flag] = eigs (AF, n ,min (10, n - 2), 'LM', opts);
eigenvalues = diag (D);
[LI2, ind] = min (eigenvalues);
extreigs.LI2 = -LI2;
inieigs.LM = V(:,ind);

function [y, normerrest, info] = newton(A, v, xi, dd, gammafoc, reltol, abstol)
%
% Compute the Newton interpolation polynomial
% nodes are in [-2,2] + d/gammafoc
maxm = size (dd, 1);
w = v;
y = w * diag (dd(1, :));
w = A (w) / gammafoc - xi(1) * w;
m = 1;
errest = w * diag (dd(2, :));
c1 = norm (errest, Inf);
y = y + errest;
w = A (w) / gammafoc - xi(m + 1) * w;
m = m + 1;
errest = w * diag (dd(m+1, :));
c2 = norm (errest(:), Inf);
while ((c1 + c2 > reltol * norm (y(:), Inf) + abstol) && ...
       (m + 1 < maxm))
  y = y + errest;
  w = A (w) / gammafoc  - xi(m + 1) * w;
% m counts the matrix-vector products
  m = m + 1;
  errest = w * diag (dd(m + 1, :));
  c1 = c2;
  c2 = norm (errest(:), Inf);
end
normerrest = c1 + c2;
if (normerrest > reltol * norm(y(:), Inf) + abstol)
% maximum number of iterations reached
  info = -maxm;
else
  info = m;
end

function [y, normerrest, info] = newtons (A, v, xi, dd, gammafoc, reltol, abstol)
%
% Compute the Newton interpolation polynomial
% nodes are in i*[-2,2] + d/gammafocfoc
%
maxm = size (dd, 1);
w = v;
y = w * diag (real (dd(1, :))); % dd(1,:) should be real
w = A (w) / gammafoc - xi(1) * w;
m = 1;
errest1 = w * diag (real (dd(2, :)));
c1 = norm (errest1(:), Inf);
wtilde = A (w) / gammafoc - real (xi(2)) * w;
errest2 = wtilde * diag (real (dd(3, :)));
c2 = norm (errest2(:), Inf);
while ((c1 + c2 > reltol * norm (y(:), Inf) + abstol) && ...
       (m + 3 < maxm))
% if the points are in -d+i*gamma*[-2,2], the divided differences are
% small. The approximation has a small initial norm, then growing. It is
% better to force some iteration. The better is the estimate on the
% number of iteration, the better is to allow a minimum number of its.
  y = y + errest1 + errest2;
  w = A (wtilde) / gammafoc - real (xi(m+1)) * wtilde + imag (xi(m+1))^2 * w;
% m counts the matrix-vector products
  m = m + 2;
  wtilde = A (w )/ gammafoc - real (xi(m + 1)) * w;
  errest1 = w * diag (real (dd(m + 1, :)));
  c1 = norm (errest1(:), Inf);
  errest2 = wtilde * diag (real (dd(m + 2, :)));
  c2 = norm (errest2(:), Inf);
end
normerrest = c1 + c2;
if (normerrest > reltol * norm (y(:), Inf) + abstol)
% maximum number of iterations reached
  info = -maxm;
else
  info = m;
end

function d = divdif(funct,k,z,h,gamma)
% z in d/gamma + [-2,2] or d/gamma + i*[-2,2]
m = length (z);
% Formula wikipedia Alternative definitions
%A = cumprod (z(:,ones(1,m)) - z(:,ones(1,m)).' + eye (m), 2);
%n = length (h);
%for i = 1:n
%  d(:,i) = funct (h(i) * gamma * z ,k); % column vector
%  d(:,i) = sum (triu (d(:,i*ones(1,m))) ./ A);
%end
%d

% Stardard recurrence
%d = funct (h * gamma * z ,k);
%for i = 2:m
%  for j = 2:i
%    d(i) = (d(i) - d(j - 1)) / (z(i) - z(j - 1));
%  end
%end

% Function matrix formulation
% COME FARE CON h diversi?
d = funct (k, z, h * gamma);

function xi = leja(m)
xi=[0.20000000000000000E+01;
 -0.20000000000000000E+01;
  0.00000000000000000E+00;
  0.11547005389072722E+01;
 -0.13174131921908798E+01;
  0.16785083483460106E+01;
 -0.17400142991341268E+01;
 -0.61122666051585473E+00;
  0.64341522576591714E+00;
  0.18859583645359135E+01;
 -0.19053465427082013E+01;
 -0.95882465993895383E+00;
  0.14252772819409354E+01;
  0.31191873036563822E+00;
 -0.15497446836687399E+01;
  0.19589552375933266E+01;
 -0.32233053684798985E+00;
 -0.19666526192556142E+01;
  0.92274120557475403E+00;
  0.17837856425002223E+01;
 -0.11437941710547146E+01;
 -0.18251194862987536E+01;
  0.12985070496776170E+01;
 -0.15963594650138674E+00;
  0.48461307997069958E+00;
 -0.14445887774072832E+01;
  0.19853372477170090E+01;
 -0.78157903185550548E+00;
  0.15642613771121558E+01;
 -0.19881134949997448E+01;
  0.79515378837945172E+00;
 -0.16527693558289136E+01;
  0.18400967618043405E+01;
 -0.47031502482404541E+00;
  0.16237742431025759E+00;
 -0.19362973057205193E+01;
  0.19283923034097246E+01;
  0.10488606867030650E+01;
 -0.10567169290471541E+01;
  0.14988429403852255E+01;
 -0.17846584511097803E+01;
 -0.12391389737003360E+01;
  0.17305853530680211E+01;
  0.39781032979286135E+00;
 -0.69722681619218396E+00;
  0.19949158007765879E+01;
 -0.18707564721468248E+01;
  0.12274622398609041E+01;
 -0.23862075208070349E+00;
 -0.16002204184230693E+01;
  0.72036969721353816E+00;
 -0.19958685914305654E+01;
  0.16242140212548106E+01;
 -0.87412702016542365E+00;
  0.19724261744339364E+01;
  0.80404158966704414E-01;
 -0.13876623788610769E+01;
  0.13631392018454864E+01;
 -0.19529498748553720E+01;
  0.98318328054576365E+00;
 -0.40340085143054333E+00;
  0.19071556630242674E+01;
 -0.16975573834583724E+01;
  0.56064317542589448E+00;
 -0.14958571038756896E+01;
  0.18124992798591713E+01;
 -0.78233177400729326E-01;
 -0.19792196201336756E+01;
  0.11048530744290983E+01;
 -0.10091469774874193E+01;
  0.19981611763672014E+01;
 -0.54661155815383133E+00;
 -0.18490785255385729E+01;
  0.15324942269853394E+01;
  0.24047293749664256E+00;
 -0.11945330474062112E+01;
  0.18638789913676468E+01;
  0.85523120854811341E+00;
 -0.19205417848388144E+01;
  0.17047215720365889E+01;
 -0.82681093074815104E+00;
 -0.16269396953652424E+01;
  0.19454358501538969E+01;
  0.44049818260554907E+00;
 -0.19985529667308248E+01;
  0.12646226928741116E+01;
 -0.13516689386557659E+01;
 -0.27959567945593367E+00;
  0.14595518634318552E+01;
 -0.17632647864745934E+01;
  0.12100991164387512E+00;
  0.19796700401009231E+01;
 -0.65578679811245655E+00;
  0.68259629139607525E+00;
 -0.18882965699713949E+01;
  0.17587184762659824E+01;
 -0.11021736407306699E+01;
  0.15977510949799085E+01;
 -0.15224243565360993E+01;
 -0.11713917398659571E+00;
  0.10156251201074136E+01;
 -0.19920609183244737E+01;
  0.19909754389942678E+01;
 -0.91656729720428121E+00;
  0.35125666439067899E+00;
  0.13342493738720118E+01;
 -0.18057573015625830E+01;
  0.19178282223700041E+01;
 -0.12783352957127079E+01;
 -0.50636888139925262E+00;
  0.88711907622909425E+00;
 -0.19728690859025828E+01;
  0.16529163699606504E+01;
  0.36324223875185128E-01;
 -0.16772535961126018E+01;
  0.18269295872608657E+01;
  0.59828071959314566E+00;
 -0.14185314526121953E+01;
  0.11892745028984089E+01;
 -0.74004480758623958E+00;
  0.19655366634970863E+01;
 -0.19448214118826375E+01;
 -0.36473230837564752E+00;
  0.13971063858841664E+01;
 -0.17201596479296610E+01;
  0.20638660864548197E+00;
  0.19993615581122488E+01;
 -0.11694111992233767E+01;
  0.76118357247152990E+00;
 -0.19840843366776040E+01;
  0.18751188672466657E+01;
 -0.15743849287764438E+01;
 -0.19837146670101830E+00;
  0.10794685548932486E+01;
 -0.18600756919982304E+01;
  0.17448849163392661E+01;
 -0.98435277450936931E+00;
  0.52115164438684380E+00;
  0.14804723469509264E+01;
 -0.58000573293431446E+00;
 -0.19994907653772023E+01;
  0.19378502518760174E+01;
 -0.14698147331485774E+01;
 -0.38665298263104081E-01;
  0.15813962348335746E+01;
 -0.19130149913849928E+01;
  0.95174241520235947E+00;
 -0.85012679882534403E+00;
  0.17978673155225322E+01;
  0.27790115008970273E+00;
 -0.12976347570633995E+01;
  0.19884036833991652E+01;
 -0.18356760596259141E+01;
 -0.43591547045755408E+00;
  0.12458760307759087E+01;
 -0.19600030758784839E+01;
  0.82352085825767629E+00;
  0.18970546328941940E+01;
 -0.10796895476808857E+01;
 -0.16398694765802584E+01;
  0.16657857586328304E+01;
  0.37556569558272346E+00;
 -0.71851159293948874E+00;
 -0.17522623869784195E+01;
  0.19528019576324605E+01;
  0.11318200384276298E+01;
 -0.13701244963593509E+01;
 -0.13843992971968599E+00;
 -0.19941800376844832E+01;
  0.13806290736436702E+01;
  0.19967087118551494E+01;
  0.62076650375226761E+00;
 -0.12179461714907480E+01;
 -0.18961260995395672E+01;
  0.59321622300415357E-01;
  0.18515979143014407E+01;
 -0.34333214126906741E+00;
  0.15174104510155204E+01;
 -0.15363309140757493E+01;
 -0.93608754429389618E+00;
  0.17169049281230784E+01;
 -0.19295036355898909E+01;
  0.74040169612773110E+00;
 -0.63351253860727552E+00;
  0.19761830519314971E+01;
 -0.17951877483432130E+01;
  0.46286329034297946E+00;
  0.13156820343408633E+01;
 -0.19974323573365131E+01;
  0.18374074686749642E+00;
  0.17720726595862875E+01;
 -0.11234112999456882E+01;
 -0.15879222569701819E+01;
  0.10318521087330463E+01;
 -0.25886694816712452E+00;
  0.19232144716501198E+01;
 -0.17087331818172968E+01;
  0.14421895057088641E+01;
 -0.80240512714215129E+00;
 -0.19760626278667983E+01;
  0.16368140642612501E+01;
  0.54134829760904868E+00;
 -0.14046511952885661E+01;
  0.19931734850658200E+01;
 -0.52540321013695324E+00;
  0.90471021449558708E+00;
 -0.18791077893082284E+01;
 -0.57496557396576657E-01;
  0.12062363247183083E+01;
 -0.10328450798122886E+01;
  0.18914741137953366E+01;
 -0.18156421615632654E+01;
  0.25994109504469465E+00;
 -0.12592408239980593E+01;
  0.15495374423760373E+01;
  0.19997771269354820E+01;
 -0.19565437514122395E+01;
  0.66510360596327245E+00;
 -0.14828859488893442E+01;
 -0.45255027071268694E+00;
  0.18197318521595538E+01;
 -0.16661085174995318E+01;
  0.11705194181332477E+01;
  0.10248896482064303E+00;
 -0.19900489631934479E+01;
  0.19687310193550804E+01;
 -0.89439995684080076E+00;
  0.16918063576731750E+01;
 -0.21660687463806830E+00;
  0.83931153119878499E+00;
 -0.17737863422802569E+01;
 -0.13342974127878295E+01;
  0.13487982074272646E+01;
 -0.67663150865144928E+00;
  0.19417324619830225E+01;
 -0.19406612086140811E+01;
  0.41871445631709769E+00;
  0.16104042052877992E+01;
 -0.16133660507757301E+01;
  0.96821857939159217E+00;
 -0.19998216287681223E+01;
  0.18577978199406866E+01;
 -0.38481592732544789E+00;
 -0.11816558948490286E+01;
  0.33020376949941965E+00;
  0.19827672639234795E+01;
 -0.76159920845964324E+00;
 -0.18429308141537173E+01;
  0.12812943279381548E+01;
 -0.17421475485852884E-01;
 -0.14561587444279338E+01;
  0.17907669681095686E+01;
  0.77832848803746946E+00;
 -0.19696673266172788E+01;
  0.14125372999016872E+01;
 -0.99696993264214240E+00;
  0.19987753701654216E+01;
 -0.17301164448261215E+01;
  0.57968332715471649E+00;
 -0.30085287689251250E+00;
  0.10660513668629648E+01;
 -0.19247788397316214E+01;
  0.19121492714876558E+01;
 -0.56471201638588608E+00;
 -0.15101018918392919E+01;
  0.17377548010153787E+01;
  0.14268736140062088E+00;
 -0.19859808012663684E+01;
  0.15079242541034690E+01;
 -0.10909202482605411E+01;
  0.19619255586965365E+01;
 -0.98568340014749664E-01;
 -0.16873779171436061E+01;
  0.11184110878802367E+01;
  0.70159825402612641E+00;
 -0.18656326928837734E+01;
 -0.12492464317924714E+01;
  0.16449731476462048E+01;
 -0.83846533209976970E+00;
  0.18342414869889936E+01;
  0.50176886337164861E+00;
 -0.15623753554246487E+01;
 -0.19967152509394768E+01;
  0.19897363639201788E+01;
 -0.48850557480153239E+00;
  0.99903538669353953E+00;
  0.22326043560717895E+00;
 -0.19007727069998448E+01;
  0.14698034001810949E+01;
 -0.14301917806537630E+01;
  0.19331056769471811E+01;
 -0.17875649221010009E+00;
 -0.94762700509045317E+00;
  0.15730360944594053E+01;
 -0.19490628039745230E+01;
  0.87169272306901247E+00;
 -0.59757074156589451E+00;
  0.18800175148995142E+01;
 -0.17899392769268649E+01;
  0.12173334522913510E+01;
 -0.13076123189439273E+01];
xi = xi(1:m);

function xi = lejas(m)
% symmetric
xi = [0.00000000000000000E+00;
 -0.20000000000000000E+01;
  0.20000000000000000E+01;
 -0.11547005389072722E+01;
  0.11547005389072722E+01;
 -0.16798869581724092E+01;
  0.16798869581724092E+01;
 -0.55212517153625673E+00;
  0.55212517153625673E+00;
 -0.18860377814061269E+01;
  0.18860377814061269E+01;
 -0.86926567048810433E+00;
  0.86926567048810433E+00;
 -0.14648353835608652E+01;
  0.14648353835608652E+01;
 -0.25358605654518529E+00;
  0.25358605654518529E+00;
 -0.19608080183933749E+01;
  0.19608080183933749E+01;
 -0.17855360051319420E+01;
  0.17855360051319420E+01;
 -0.13076036318707878E+01;
  0.13076036318707878E+01;
 -0.70927408113420620E+00;
  0.70927408113420620E+00;
 -0.19857559840157573E+01;
  0.19857559840157573E+01;
 -0.15746119124852025E+01;
  0.15746119124852025E+01;
 -0.38641404790579925E+00;
  0.38641404790579925E+00;
 -0.10252007303405923E+01;
  0.10252007303405923E+01;
 -0.19238218904644926E+01;
  0.19238218904644926E+01;
 -0.11356505078523751E+00;
  0.11356505078523751E+00;
 -0.18328648391574527E+01;
  0.18328648391574527E+01;
 -0.13869875745457783E+01;
  0.13869875745457783E+01;
 -0.17299132735089220E+01;
  0.17299132735089220E+01;
 -0.79110322966313817E+00;
  0.79110322966313817E+00;
 -0.19950094737219035E+01;
  0.19950094737219035E+01;
 -0.12255720758043205E+01;
  0.12255720758043205E+01;
 -0.46869711976191569E+00;
  0.46869711976191569E+00;
 -0.16249179463353562E+01;
  0.16249179463353562E+01;
 -0.19439167904762229E+01;
  0.19439167904762229E+01;
 -0.95416908607531581E+00;
  0.95416908607531581E+00;
 -0.18233461816615804E+00;
  0.18233461816615804E+00;
 -0.18604221037194564E+01;
  0.18604221037194564E+01;
 -0.15165949191908132E+01;
  0.15165949191908132E+01;
 -0.63208883307573149E+00;
  0.63208883307573149E+00;
 -0.19757001452981686E+01;
  0.19757001452981686E+01;
 -0.10945725314996602E+01;
  0.10945725314996602E+01;
 -0.17583097219194987E+01;
  0.17583097219194987E+01;
 -0.32085293148252758E+00;
  0.32085293148252758E+00;
 -0.13487547756021807E+01;
  0.13487547756021807E+01;
 -0.19982428230270117E+01;
  0.19982428230270117E+01;
 -0.51257930118453976E-01;
  0.51257930118453976E-01;
 -0.19054515139305148E+01;
  0.19054515139305148E+01;
 -0.14300340970675731E+01;
  0.14300340970675731E+01;
 -0.91121106035049493E+00;
  0.91121106035049493E+00;
 -0.16533046750271887E+01;
  0.16533046750271887E+01;
 -0.59277547325880930E+00;
  0.59277547325880930E+00;
 -0.18114456322046653E+01;
  0.18114456322046653E+01;
 -0.12609004281591163E+01;
  0.12609004281591163E+01;
 -0.19904989930840369E+01;
  0.19904989930840369E+01;
 -0.75108525679913374E+00;
  0.75108525679913374E+00;
 -0.15471104012629255E+01;
  0.15471104012629255E+01;
 -0.19526142499017065E+01;
  0.19526142499017065E+01;
 -0.42646412321259564E+00;
  0.42646412321259564E+00;
 -0.10611747892375527E+01;
  0.10611747892375527E+01;
 -0.17068490728341663E+01;
  0.17068490728341663E+01;
 -0.21711234979315699E+00;
  0.21711234979315699E+00;
 -0.18737209418635996E+01;
  0.18737209418635996E+01;
 -0.11885139040381980E+01;
  0.11885139040381980E+01;
 -0.19695392132630480E+01;
  0.19695392132630480E+01;
 -0.83065749929225863E+00;
  0.83065749929225863E+00;
 -0.15998430128521397E+01;
  0.15998430128521397E+01;
 -0.50969282235341318E+00;
  0.50969282235341318E+00;
 -0.19329922003399136E+01;
  0.19329922003399136E+01;
 -0.14887639818661527E+01;
  0.14887639818661527E+01;
 -0.82500260421577265E-01;
  0.82500260421577265E-01;
 -0.19993753417027305E+01;
  0.19993753417027305E+01;
 -0.99031285197543495E+00;
  0.99031285197543495E+00;
 -0.17984694537070336E+01;
  0.17984694537070336E+01;
 -0.12856068973614110E+01;
  0.12856068973614110E+01;
 -0.67083107410985932E+00;
  0.67083107410985932E+00;
 -0.18467090979098812E+01;
  0.18467090979098812E+01;
 -0.29086225313969671E+00;
  0.29086225313969671E+00;
 -0.19813115154694114E+01;
  0.19813115154694114E+01;
 -0.14082128486138015E+01;
  0.14082128486138015E+01;
 -0.11248472454724956E+01;
  0.11248472454724956E+01;
 -0.17439345363407788E+01;
  0.17439345363407788E+01;
 -0.35608248193825248E+00;
  0.35608248193825248E+00;
 -0.19141794488578103E+01;
  0.19141794488578103E+01;
 -0.16393881747858858E+01;
  0.16393881747858858E+01;
 -0.89037695394600458E+00;
  0.89037695394600458E+00;
 -0.19930460273640178E+01;
  0.19930460273640178E+01;
 -0.14722764520692921E+00;
  0.14722764520692921E+00;
 -0.13303270238321043E+01;
  0.13303270238321043E+01;
 -0.18954957715477345E+01;
  0.18954957715477345E+01;
 -0.73002695721571742E+00;
  0.73002695721571742E+00;
 -0.15322007488124774E+01;
  0.15322007488124774E+01;
 -0.53093030417538989E+00;
  0.53093030417538989E+00;
 -0.19652740490686988E+01;
  0.19652740490686988E+01;
 -0.16937068017621355E+01;
  0.16937068017621355E+01;
 -0.10432560583884489E+01;
  0.10432560583884489E+01;
 -0.23104523884455920E-01;
  0.23104523884455920E-01;
 -0.19969351840405918E+01;
  0.19969351840405918E+01;
 -0.12072030165460474E+01;
  0.12072030165460474E+01;
 -0.17725812834896795E+01;
  0.17725812834896795E+01;
 -0.14479963855736195E+01;
  0.14479963855736195E+01;
 -0.44735129673309959E+00;
  0.44735129673309959E+00;
 -0.18233260155856839E+01;
  0.18233260155856839E+01;
 -0.81116207966067522E+00;
  0.81116207966067522E+00;
 -0.19387483888387484E+01;
  0.19387483888387484E+01;
 -0.97174274699549334E+00;
  0.97174274699549334E+00;
 -0.15872920982131618E+01;
  0.15872920982131618E+01;
 -0.61287752274714813E+00;
  0.61287752274714813E+00;
 -0.19880574665310511E+01;
  0.19880574665310511E+01;
 -0.13681257452441069E+01;
  0.13681257452441069E+01;
 -0.23539334410327134E+00;
  0.23539334410327134E+00;
 -0.18671540015844768E+01;
  0.18671540015844768E+01;
 -0.11699972963064553E+01;
  0.11699972963064553E+01;
 -0.16671457859178513E+01;
  0.16671457859178513E+01;
 -0.19997819819489917E+01;
  0.19997819819489917E+01;
 -0.33879692067760181E+00;
  0.33879692067760181E+00;
 -0.93193510956905001E+00;
  0.93193510956905001E+00;
 -0.17193556727601473E+01;
  0.17193556727601473E+01;
 -0.19564712782337348E+01;
  0.19564712782337348E+01;
 -0.68810730108543350E+00;
  0.68810730108543350E+00;
 -0.15017426944300611E+01;
  0.15017426944300611E+01;
 -0.13044591694027682E+00;
  0.13044591694027682E+00;
 -0.19189501051074362E+01;
  0.19189501051074362E+01;
 -0.12449462828393627E+01;
  0.12449462828393627E+01;
 -0.18400500193660232E+01;
  0.18400500193660232E+01;
 -0.48936152483709394E+00;
  0.48936152483709394E+00;
 -0.11091391592482718E+01;
  0.11091391592482718E+01;
 -0.19785744492297805E+01;
  0.19785744492297805E+01;
 -0.15609602234431856E+01;
  0.15609602234431856E+01;
 -0.77167718033496546E+00;
  0.77167718033496546E+00;
 -0.17656073238286094E+01;
  0.17656073238286094E+01;
 -0.66682609227377843E-01;
  0.66682609227377843E-01;
 -0.14186948906067163E+01;
  0.14186948906067163E+01;
 -0.19960372860004174E+01;
  0.19960372860004174E+01;
 -0.57293757345009144E+00;
  0.57293757345009144E+00;
 -0.18807232956436515E+01;
  0.18807232956436515E+01;
 -0.10084464422900927E+01;
  0.10084464422900927E+01;
 -0.16133381879660718E+01;
  0.16133381879660718E+01;
 -0.40436849775440142E+00;
  0.40436849775440142E+00;
 -0.19481910481270022E+01;
  0.19481910481270022E+01;
 -0.12967412655742256E+01;
  0.12967412655742256E+01;
 -0.85067147549021660E+00;
  0.85067147549021660E+00;
 -0.18048858729861337E+01;
  0.18048858729861337E+01;
 -0.27280542126864110E+00;
  0.27280542126864110E+00;
 -0.19727206400250756E+01;
  0.19727206400250756E+01;
 -0.14763121584831804E+01;
  0.14763121584831804E+01;
 -0.10782667121173539E+01;
  0.10782667121173539E+01;
 -0.19006344714612742E+01;
  0.19006344714612742E+01;
 -0.16697056908903357E+00;
  0.16697056908903357E+00;
 -0.16869141167379311E+01;
  0.16869141167379311E+01;
 -0.65170298750894351E+00;
  0.65170298750894351E+00;
 -0.19988653296263976E+01;
  0.19988653296263976E+01;
 -0.13583663084152575E+01;
  0.13583663084152575E+01;
 -0.11407374157356993E+01;
  0.11407374157356993E+01;
 -0.17373227019455160E+01;
  0.17373227019455160E+01;
 -0.19837088666911726E+01;
  0.19837088666911726E+01;
 -0.37147733137072370E+00;
  0.37147733137072370E+00;
 -0.12721419615906764E+01;
  0.12721419615906764E+01];
xi = xi(1:m);

function d = phiktaylordd (k, xi, gamma, s)
% d = phiktaylordd (k, xi, gamma)
%
% Compute accurate divided differences for the phi_k functions at nodes
% xi(:). If xi(:) are in [-2,2]+c/gamma, then they are the divided 
% differences for phi_k (c + gamma * xi) at nodes xi(:). All the
% divided differences for phi_0, ..., phi_k are computed and put in the 
% columns of d. That is, d corresponds to
%
% m = length (xi);
% for l = 0:k
%   for i = 1:m
%     d(i, l + 1) = phik (l, gamma * xi(i));
%     for j = 1:i - 1
%       d(i, l + 1) = (d(i, l + 1) - d(j, l + 1)) / (xi(i) - xi(j));
%     end
%   end
% end
%
l = max (gamma * (abs (xi(1:3)) + 1));
if (nargin == 3)
  s = floor (2 * l);
end
m = length (xi);
gamma = gamma / s;
dg = zeros(m, k);
Fk = zeros(m + k, m + k);
t = toeplitz (0:m - 1);
F = gamma .^ t ./ factorial (t);
A = (0:m - 1)';
B = 1:k;
g = gamma .^ A(:, ones(1, k)) ./ ...
    factorial (A(:, ones(1, k)) + B(ones(m, 1), :));
d = F(:,1);
dg = g;
error = 1;
l = 1;
xi = xi(:).';
% pointers to main diagonal
iddiag = (1:m + 1:m ^ 2);
%idtril = (2:m);
%idtriu = (m+1:m:m^2-m+1);
%for i = 2:m-1
%  idtril = [idtril,(i-1)*m+1+i:i*m];
%  idtriu = [idtriu,i*m+i:m:m^2-m+i];
%end
while (error > 0)
  F(iddiag) = (gamma / l) * F(iddiag) .* xi;
  g(1, :) = gamma * xi(1) * g(1, :) ./ (l + (1:k));
  idx1 = iddiag;
  for i = 1:m-1
% pointers to sup diagonals
    idx2 = idx1(1:m - i);
    idx1 = idx1(2:m - i + 1) - 1;
    F(idx1) = (gamma / (l + i)) * (F(idx1) .* xi(i + 1:m) + F(idx2));
    g(i + 1, :) = gamma * (xi(i + 1) * g(i + 1, :) + g(i, :)) ./ ...
	(l + i + (1:k));
  end  
  F = F + triu (F).' - diag (F(iddiag));
%  F(idtril) = F(idtril) + F(idtriu);
  dg = dg + g;
  error = max (abs (F(2:m, 1) - d(2:m)) ./ abs (F(2:m, 1)));
  d = F(:, 1);
  l = l + 1;
end
F(iddiag) = exp (gamma * xi);
Fk(1:m, 1:m) = tril (F);
if (k > 0) 
% there is a bug in Matlab, giving error for toeplitz(zeros([0:-1]),
% instead of empty matrix
  t = toeplitz (0:k - 1);
  Fk(m + 1:m + k, m + 1:m + k) = triu (ones(k) ./ factorial (t) ./ s .^ t);
% compute the phik functions in gamma * xi(1)
  if (abs (gamma * xi(1)) > 0.5)
    dg(1,1) = (exp (gamma * xi(1)) - 1) / (gamma * xi(1));
    for i = 2:k
      dg(1, i) = (dg(1, i - 1) - 1 / factorial (i - 1)) / (gamma * xi(1));
    end
  else
    c = 1 ./ factorial (1:k + 14);
    for i = 1:k
      dg(1, i) = polyval (fliplr (c(i:i + 14)), gamma * xi(1));
    end			  
  end
end
Fk(1:m, m + 1:m + k) = dg / diag (s .^ (1:k));
%Fk(1:m, m+1:m + k) = bsxfun (@rdivide, dg, s .^ (1:k));
d = Fk(:, [1, m + 1:m + k]);
for i = 1:s-1
  d = Fk * d;
end
d = d(1:m,:);

%!test
%! m = 10;
%! xi = leja (10);
%! c = -4;
%! gamma = 3;
%! xi = xi + c / gamma;
%! k = 2;
%! for l = 0:k
%!   for i = 1:m
%!     d(i, l + 1) = phik (l, gamma * xi(i));
%!     for j = 1:i - 1
%!       d(i, l + 1) = (d(i, l + 1) - d(j, l + 1)) / (xi(i) - xi(j));
%!     end
%!   end
%! end
%! d1 = phiktaylordd (k, xi, gamma);
%! assert (max (abs (d - d1)) < 1e-14)