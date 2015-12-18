function [w, stats] = expmvp(t, A, u, tol, m)
% EXPMVP - Evaluates a linear combinaton of the phi functions
%          evaluated at tA acting on vectors from u, that is 
%
%          w = phi_0(tA) u(:, 1) + t phi_1(tA) u(:, 2) + 
%              t^2 phi_2(tA) u(:, 3) + ...  
%
%          The evaluation expresses eveything in terms of the highest
%          order phi function and evaluates the action of this on a
%          vector using a Krylov technique and then computes w using
%          the recurrence relation.
%
%          The size of the Krylov subspace is changed dynamically
%          during the integration. The Krylov subspace is computed
%          using Arnoldi if A is non-symmetric and Lancozs if A is
%          symmetric. 
%
% PARAMETERS:
%   t    - constant value represent time.
%   A    - the matrix argument of the phi functions.
%   u    - the matrix with columns representing the vectors to be
%          multiplied by the phi functions.
%   tol  - the convergence tolarance required. 
%   m    - an estimate of the appropriate Krylov size.
%
% RETURNS:
%   w        - the linear combination of the phi functions
%              evaluated at tA acting on the vectors from u.
%   stats(1) - number of substeps
%   stats(2) - number of rejected steps
%   stats(3) - number of Krylov steps
%   stats(4) - number of matrix exponentials

% n is the size of the original problem
% p is the number of phi functions

% NOTE: some blind modifications has been made by Seb. This version it 
% differs from the original one found online at 
% http://www1.maths.leeds.ac.uk/~jitse/expmvp.m by Jitse Niesen
% See OLD commented and NEW uncommented.

[n, ppo] = size(u);
p = ppo - 1;
nnze = nnz(A);

% Compute scaling factor
if ppo > 2 nu = 2^(-ceil(log2(norm(u, 1)))); else nu = 1; end

% Initial condition
w = u(:, 1);

% Flip the rest of the u matrix
u = nu*fliplr(u(:, 2:end));

% construct the augmented matrix A
if p == 1; K = 0; else K = spdiags(ones(p, 1), 1, p, p); end
A = [A, u; zeros(p, n), K];

% Check inputs
if nargin < 5
  m = 10;
  if nargin < 4
    tol = 1.0e-7;
    if nargin < 3
      error('phipm:NotEnoughInputs',...
            'Not enough input arguments.');
    end  
  end
end

% Krylov parameters
mmax = 100;
m_new = m;

% Preallocate matrix
V = zeros(n + p, mmax + 1); 

% Initializing the variables
step = 0; 
krystep = 0;
ireject = 0;
reject = 0;
exps = 0;
happy = 0;
sgn = sign(t);
t_now = 0; 
t_out = abs(t);
j = 0;

% Compute and initial starting approximation for the timestep
tau = t_out;

% Setting the safety factors and tolarance requirements
gamma = 0.6; 
delta = 1.4; 

% Used in the adaptive selection
oldm = NaN; oldtau = NaN; omega = NaN;
orderold = true; kestold = true;

% Iterate until we reach the final time t
while t_now < t_out

  % Create initial w
  tt = flipud(cumprod(t_now./(1:p-1)'));
  w = [w; tt/nu; 1/nu];

  % Compute necessary starting information
  if j == 0
  
    % Initialize the matrices V and H
    H = zeros(mmax + 1 + p, mmax + 1 + p);

    % Normalize initial vector
    beta = norm(w);
    if beta == 0
      
      % Multiplying with a zero vector, hence result is zero
      % Finish all in one step
      reject = reject + ireject;
      step = step + 1;
      tau = t_out - t_now;
      break;
    
    end;

    % The first Krylov basis vector
    
    % OLD
     %V(:, 1) = w ./ beta;
    
    % NEW modified by Seb
    V(:, 1) = w(1:end-1) ./ beta;
  
  end;

  % Not symmetric use the Arnoldi process
  while j < m

    % Arnoldi
    j = j + 1;
    vv = A * V(:, j);
     for i = 1:j
      H(i, j) = V(:, i)' * vv;
      vv = vv - H(i, j) * V(:, i);
    end
    krystep = krystep + 1;
    s = norm(vv); 
    
    % Happy breakdown
    if s < tol
      happy = 1;
      tau = t_out - t_now;
      break;
    end;
    H(j + 1, j) = s;
    V(:, j + 1) = vv ./ s;
    
  end   
    
  % Keep a record of H
  H2 = H;
 
  % We use the vector e1 in the computations
  H(1, j + 1) = 1; 
  
  % Construct the augmented matrix
  h = H(j + 1, j); 
  H(j + 1, j) = 0;
  
  % Compute the exponential of the augmented matrix
  
  % OLD
  %[F,hnorm] = expmnorm(sgn*tau*H(1:j+p, 1:j+p));  
  
  %NEW - seb modified
  F = expm(sgn*tau*H(1:j+p, 1:j+p)); 
  hnorm = norm(F);
  
  exps = exps + 1;
  
  % Local truncation error estimation
  % OLD
  %err = abs(beta * h * F(j, j + 1));
  
  % new modified by Seb
  err = abs(beta * h * F(j, j));
  
  % Error per unit step
  oldomega = omega;
  omega = t_out * err / (tau * tol);

  % Estimate order
  if m == oldm && tau ~= oldtau && ireject >= 1
    order = max(1, log(omega/oldomega) / log(tau/oldtau));
    orderold = false;
  elseif orderold || ireject == 0
    orderold = true;
    order = j/4;
  else
    orderold = true;
  end;
  % Estimate k
  if m ~= oldm && tau == oldtau && ireject >= 1
    kest = max(1.1, (omega/oldomega) ^ (1/(oldm-m)));
    kestold = false;
  elseif kestold || ireject == 0
    kestold = true;
    kest = 2;
  else
    kestold = true;
  end;
  
  % This if statement is the main difference between fixed and
  % variable m  
  oldtau = tau; oldm = m;
  if happy == 1

    % Happy breakdown; wrap up
    omega = 0;
    tau_new = tau;
    m_new = m;

  elseif j == mmax && omega > delta
  
    % Krylov subspace to small and stepsize to large
    tau_new = tau * (omega / gamma) ^ (-1 / order); 

  else
  
    % Determine optimal tau and m
    tau_opt = tau * (omega / gamma) ^ (-1 / order);
    m_opt = max(1, ceil(j + log(omega / gamma) / log(kest)));
    nom = 5 + max(log(hnorm), 0) / log(2); % number of mult's in expm
    
    % Cost of Arnoldi
    cost1 = (j * nnze + (j^2 + 5) * (n + p) + nom * j^3) * ...
            ceil((t_out-t_now) / tau_opt); 
    cost2 = (m_opt * nnze + (m_opt^2 + 5) * (n + p) + nom * m_opt^3) ...
            * ceil((t_out-t_now) / tau); 
   
    % Determine whether to vary tau or m
    if cost1 < cost2 
      tau_new = tau_opt;
      m_new = m;
    else
      m_new = m_opt;
      tau_new = tau;
    end;
  
  end;
  
  % Check error against target
  if omega <= delta
     
    % Yep, got the required tolerance; update 
    reject = reject + ireject;
    step = step + 1;
    
    % Using the corrected quantity
    
    % OLD
    %F(j + 1, j) = h * F(j, j + 1);
    %w = beta * V(1:n, 1:j + 1) * F(1:j + 1, 1);
    
    % NEW:
    F(j, j) = h * F(j, j);
    w = beta * V(1:n, 1:j) * F(1:j, 1);
    
    
    % Update t
    t_now = t_now + tau;
    j = 0;
    ireject = 0;
  
  else
  
    % Nope, try again
    H = H2;
    ireject = ireject + 1; 
  
  end;
  
  % Safety factors for tau
  tau = min(t_out - t_now, max(tau/5, min(5*tau, tau_new)));
    
  % Safety factors for m
  m = max(1, min(mmax, max(floor(3/4*m), min(m_new, ceil(4/3*m)))));
  
end
stats = [step, reject, krystep, exps]; 