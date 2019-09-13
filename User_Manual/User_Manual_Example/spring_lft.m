function [lftfun] = spring_lft(k1,m,k2lim,Clim)
%definition of the number of inputs, states, outputs and omega of the LFT
u_count = 1;
x_count = 2;
y_count = 1;
om_count = 1;

%definition of the matrix delta 
DeltaSym = {
    'parName' 'indStartDiag' 'indStopDiag' 'LowerBound' 'HigherBound' 'toIdentify', 'lb', 'ub'; 
    'c' 1 1 Clim(1) Clim(2) 1 -1 1;
    'k2' 2 2 k2lim(1) k2lim(2) 1 -1 1;
        };
delta_count = 0;
for i=2:size(DeltaSym,1)
    delta_count = delta_count + DeltaSym{i,3} - DeltaSym{i,2} + 1;
end

%definition of the matrix theta (= ZETA)
if om_count == 1
    om = sym(['om1_1']);
else
om = sym('om', om_count);
end

theta_sym = [
    om(1)*om(1)*om(1);
    ];
theta_count = size(theta_sym,1);

%construction of the matrices of the LTI part
LTI = struct(...
    'A', zeros(x_count,x_count),...
    'B1', zeros(x_count,delta_count),...
    'B2', zeros(x_count,theta_count),...
    'B3', zeros(x_count,u_count),...
    'C1', zeros(delta_count,x_count),...
    'D11', zeros(delta_count,delta_count),...
    'D12', zeros(delta_count,theta_count),...
    'D13', zeros(delta_count,u_count),...
    'C2', zeros(om_count,x_count),...
    'D21', zeros(om_count,delta_count),...
    'D22', zeros(om_count,theta_count),...
    'D23', zeros(om_count,u_count),...
    'C3', zeros(y_count,x_count),...
    'D31', zeros(y_count,delta_count),...
    'D32', zeros(y_count,theta_count),...
    'D33', zeros(y_count,u_count));

%definition of the matrices
% dx1 = x2
LTI.A(1,2) = 1;
% dx2 = -k1/m*x1-1/m*w1-1/m*w2+1/m*u
LTI.A(2,1) = -k1/m;
LTI.B1(2,1) = -1/m;
LTI.B1(2,2) = -1/m;
LTI.B3(2,1) = 1/m;

% z1 = x2
LTI.D12(2,1) = 1;
% z2 = zeta
LTI.C1(1,2) = 1;

% om1 = x1
LTI.C2(1,1) = 1;
% LTI.C2(2,1) = 3;

% y1 = x1
LTI.C3(1,1) = 1;

%initial definition of the LFTfun
lftfun = struct(...
    'LTI', LTI,...
    'DeltaSym', {DeltaSym},...
    'DeltaVal', zeros(delta_count,delta_count)...
        );
    
% save additional data for reference purpose
lftfun.theta_sym = theta_sym;
lftfun.u_count = u_count;
lftfun.x_count = x_count;
lftfun.y_count = y_count;
lftfun.om_count = om_count;
lftfun.theta_count = theta_count;
lftfun.delta_count = delta_count;

%function to finalize the definition of LFTfun
lftfun = lft_finalize(lftfun);

end