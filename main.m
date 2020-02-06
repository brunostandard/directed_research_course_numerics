%% -------------------------Parameters and setup-------------------------
clear
format long
N       = 81;  % number of space steps
x0      = 0;    % start (space)
xf      = 5;    % end (space)
t0      = 0;    % start (time)
tf      = 20;   % end (time)
M       = 161;   % number of time steps
dx      = linspace(x0,xf,N); % mesh points over space
dt      = linspace(t0,tf,M); % mesh points over time
h 		  = (dx(2)-dx(1));  % length of steps in space
k 		  = (dt(2)-dt(1));  % length of time step
Ja      = -1; % flux
c       = @(t) 1; %r clean up rate == kappa
d       = @(t) 1;   % our rate of oligomer formation == nu
%%
% ------------------ sovling for first PDE --------------------
% A_t = A_{xx} - c(t) A - d(t) A ^2 
% with the flux condtion A_x(0,t) = -Ja

%
a       = zeros(N,M);   % rows = space ; columns = time
a       = evolve_diff(a,dx,dt,c,d,Ja); % solve
% next solve the non-linear part: a_t = f(a(x,t))
%evolve_src(a,f,dx,dt);
hold on;

%% --------------------------- Stead state solution plots ----------------------
% Next calculate with steady-state solution for the following equation
% A_xx = c*A + d*A^2. Assume that c and d are constant. 
% This equation can be reduced to the following:
% A_x = -(c(t)*A^2 + d(t)*(2/3)*a^3)^(1/2).
% Let A(x) = a b(z) where z = sqrt(c)*x and a = c/d so that the above can be re-written
% as   b_xx = -b - b^2 and hence b_x = - (b^2 + (2/3)*b^3)^(1/2)

%using the flux condition a_x|(x=0) = -1;
T 		= 200;    % the number points in space, for plotting
x 		= linspace(0,xf,T); % our mesh points in space
c1 		= 0.1071475644353; % the number corresponding to the flux condition, derived
                      % numerically
% our analytical solution:
b 		= @(x) -3/2 + (3/2)*((1+c1*exp(-x))/(1-c1*exp(-x) ))^2;  
A     = @(x) c(1)*b(x/sqrt(c(1)))/d(1);

fprintf('A"(0) = %g \n', (A(0.01)-A(0))/0.01);
fprintf('a"(0) = %g \n', (a(2,end)-a(1,end))/h);
%produce a vector of y or "A" values for plot.
for i = 1:T
  y(i) 	= A(x(i));
end
% plot the steady state solution
w3 		= plot(x,y,'k', 'linewidth', 3.5);
legend(w3,'steady-state');


sum_diff    = 0;
for i = 1:N
  diff(i) = abs(a(i,end)-A(dx(i)));
  sum_diff += diff(i);
end
avg_diff = sum_diff/N;
fprintf('h = %g \nk = %g \nsum difference = %g \n', h,k, avg_diff);


%%
% ---------------------Solution for just the heat equation with flux-------
%u_t = u_{xx} with u(x,0) = 0 and u_x(0,t) = -j

u     = @(x,t)  -x*erf(x/(2*sqrt(t))) - (2*sqrt(t))*exp(-x^2/(4*t))/sqrt(pi) + x;
t     = linspace(0.01,tf,M);
x     = linspace(0.01,xf,N);
figure('name', 'Heat Eqn w/ Flux condition');
for i = 1:M
  for j = 1:N
    y2(j,i) = Ja*u(x(j),t(i));
  end
  plot(x,y2(:,i),'k');
  hold on;
end
title('u_t = u_{xx}  , u_x(0,t) = -j, & u(x,0) = 0');
xlabel('x [space]');
ylabel('u [concentration]');



%%
%----------------------now for the other PDE --------------
%G_t = G_xx + (1/2) A^2 - G
%We will work with the linear part separately from the non-linear part. 
% G_t = G_xx first with 

g       = zeros(N,M);   % the matrix solution (prep)
%g(:,1)  = 0.1*ones(N,1);  % assume we have some bad AB to begin with

e       = @(t) 1; % our clearance rate == mu

g       = evolve_diff_G(g,a,dx,dt,d,e); % solving for linear part;
%%
% ------------------------Steady state for other equation
%
figure('name', 'steady-state Bad AB');
dydt    = @(x,y) [y(2); y(1)-(1/2)*A(x)^2];
[x,y]   = ode45(dydt, dx, [0;0]);

w1 = plot(y(1,1),y(1,2),'o-b','linewidth',3);
hold on;
w2 = plot(y(end,1),y(end,2),'o-k','linewidth',3);
legend([w1;w2], 't = 0', 't = tf');
plot(y(:,1),y(:,2),'o-r');
hold on;
title('Steady state solution for Bad AB equation');
xlabel('g_1 == G');
ylabel('g_2 == G_{x}');


%----------------------now solving for health ODE ------------------------
%Hdot(t) = - o_n H * G(0,t);	// if apoptosis
%Hdot(t) = - o_v H^(2/3) * G(0,t); // if volume
% the solutions can be done analytically, 
%H = exp(-o_n*G(0,t)*t) for the first version of the ode, and
%
%H = (2*o_v*G(0,t)+1)^(-1) for the second version of the ode. 
% assume:
o_n     = 1;
o_v     = 1;

H_1       = @(i,t) exp(-o_n*g(1,i)*dt(i)); % exact solution for first version
H_2       = @(i,t) 1/(2*o_v*g(1,i)*t+1);  % exact solution for second version
w       = zeros(M,2);
for i = 1:M
  w(i,1) = H_1(i,dt(i));
  w(i,2) = H_2(i,dt(i));
end
% plot the steady state solution
figure('name','health vs time');
w1 = plot(dt,w(:,1),'k', 'linewidth', 3.5);
hold on;
w2 = plot(dt,w(:,2),'b', 'linewidth', 3.5);
legend([w1;w2],'if apoptosis','if volume');
title('Health vs time');
xlabel('t [time]');
ylabel('H [health]')

