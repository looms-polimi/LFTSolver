clear all, close all

%definition of time
% t=linspace(0,100,101);

%definition of the input
% u=1*t;%*sin(t*5);


%definition of the initial state
x0 = [0 0 0];

%computation
sol = ode15s(@(t,x)m_k_c_nonlinear(t,x),[0 70],x0);
t=sol.x';x1=sol.y(1,:)';x2=sol.y(2,:)'; u=sol.y(3,:)';

figure;
plot(t,x1)
figure;
plot(t,x2)
figure;
plot(t,u)

F=[t u]
lft_y=x1;