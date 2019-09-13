function dxdt = m_k_c_nonlinear (t,x)

%definition of the parameters
b=8; %N*sec/m
k1=20; %N/m
k2=13; %N/m
m=1; %kg
K=500;

dxdt = zeros(3,1);
dxdt(1)=x(2);
dxdt(2)=-(k1*x(1))/m-(k2*x(1)^3)/m-(b*x(2))/m+x(3)/m;
dxdt(3)=K*0.1*cos(0.1*t);

end