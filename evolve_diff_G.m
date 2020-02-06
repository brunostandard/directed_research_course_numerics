function g = evolve_diff_G(g,a,dx,dt,d,e)
N         = size(dx,2);   % number of space steps
M         = size(dt,2);   % number of time steps
h         = dx(2)-dx(1);  % length of space step
k         = dt(2)-dt(1);  % length of time step

%///////////////////////////////////////////////////////////
%         	Constructing L-Matrix and I-matrix
%///////////////////////////////////////////////////////////
v       = ones(N,1);    % helper vector
B       = [v, -2*v, v];   % helper matrix
r       = [-1 0 1];       % ordering for spdiagonal matrix
L       = spdiags(B,r,N,N);  % tri-diagonal matrix
L(1,2)  = 2;

I       = eye(N); 		% just identity matrix

% our iterative method is 
% A(J+1) = (I - (k/h^2) * L )^-1 (A(j)) 

%///////////////////////////////////////////////////////////
%         	Generating solutions
%///////////////////////////////////////////////////////////

%
%f       = @(t,G,A) + d(t)*A.*A - e(t)*G; % the non-linear part of the PDE
%figure('name','normal AB vs space');
%for j = 1:M-2
%  gbar
%  g(:,j+1) = (I - (k/h^2)*L)\g(:,j);
%end


figure('name','normal AB vs space');
for j = 1:M-1
  gbar      = (I - (k/h^2)*L)\g(:,j);
  g0    = g(:,j);
  f     = @(g) k*(-e(dt(j))*g  + (1/2)*d(dt(j))*a(:,j).*a(:,j))+gbar;
  df    = @(g) k*(-e(dt(j)));
  g(:,j+1) = newton_method(f,df,g0,k);
  g(:,j+1) = f(g(:,j+1));
  plot(dx,g(:,j),'b','linewidth',0.5);
  %pause(0.1);
  hold on;
end
title('G_t = G_{xx} + (1/2) d(t) A^2 - e(t) G');
xlabel('x [space]');
ylabel('G [BAD AB concentration]');