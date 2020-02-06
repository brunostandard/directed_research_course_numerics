function a = evolve_diff(a,dx,dt,c,d,Ja)
N         = size(dx,2);   % number of space steps
M         = size(dt,2);   % number of time steps
h         = dx(2)-dx(1);  % length of space step
k         = dt(2)-dt(1);  % length of time step

%///////////////////////////////////////////////////////////
%         	Constructing L-Matrix and I-matrix
%///////////////////////////////////////////////////////////
v       = ones(N,1);    % helper vector
B       = [v, -2*v, v];   % helper matrix
e       = [-1 0 1];       % ordering for spdiagonal matrix
L       = spdiags(B,e,N,N);  % tri-diagonal matrix
L(1,2)  = 2;

%L

I       = eye(N); 		% just identity matrix

%Ja      = -1;            % the (absolute value of the) flux
V       = [Ja; zeros(N-1,1)]; % auxillary vector

% our iterative method is 
% A(J+1) = (I - (k/h^2) * L )^-1 (A(j) - 2*k/h * V) 

%///////////////////////////////////////////////////////////
%         	Generating solution
%///////////////////////////////////////////////////////////

figure('name','normal AB vs space');
for j = 1:M-1
  abar      = (I - (k/h^2)*L)\(a(:,j)- 2*(k/h)*V);
  a0    = a(:,j);
  f     = @(a) k*(-c(dt(j))*a  - d(dt(j))*a.*a)+abar;
  df = @(a,t) k*(-c(dt(j))-2*a*d(dt(j)));
  a(:,j+1) = newton_method(f,df,a0,k);
  a(:,j+1) = f(a(:,j+1));
  plot(dx,a(:,j),'r','linewidth',0.5);
  %pause(0.1);
  hold on;
end

title('A_t = A_{xx} - c(t) A - d(t) A^2');
xlabel('x [space]');
ylabel('A [AB concentration]');

end