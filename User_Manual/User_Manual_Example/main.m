
clc, clear all, close all

%import input (Force 'F') NB every input or state comes with its timevector
load F.mat;
t = F(:,1); 
F = F(:,2);

%import the output (position 'x') 
load('lft_y');

%definition of the constant and known parameters
k1=20; %[kN/m]
m=1; %[Kg]

%definition of the limits of the uncertain parameters
%in this case the values to be found is 'k2=13' and 'c=8'
k2lim = [1 200];
Clim = [1 200];

%load lft function and solver options
lftfun = spring_lft(k1,m,k2lim,Clim);

%load LFT solver options
LFTsolverOptions = lftSet('RelTol', 1e-4,...
    'AbsTol', 1e-8,...
    'SolutionInterpMethod', 'spline',...
    'SolutionTimeSpan', t,...
    'SensAlgorithm', 'ode15s',...
    'OversamplingMethod', 'spline' );

%define the input
Input = struct('Type', 'interpolated', 'Samples', F, 'Time', t);

%define the initial conditions of the model's state variables
InitialConditions = struct('StateInitialConditions', [0; 0]);

%load LFT estimator options
LFToptimOptions = lftOptSet('Display', 'iter',...
    'MaxIter', 120, ...
    'TolFun', 1e-8 ...
    );

% estimate
% initial guess of parameters in normalized form
lftfun.DeltaVal = diag([0 0]);

%estimation
[DELTA_opt,fval,J,grad,H,CN, history] = lftOptDelta(lftfun,Input,InitialConditions,LFTsolverOptions,lft_y,LFToptimOptions);

norm2abs(lftfun, lftfun.DeltaVal, DELTA_opt);

H=H(:,:,end), CN=CN(end)



