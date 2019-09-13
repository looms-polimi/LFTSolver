function [Input fval CostFun gradient ConditionNumber history, varargout] = lftOptInput(varargin)



    %----------------------------------------------------------------------
    % Global Variables
    %----------------------------------------------------------------------
    global J;
    global G;
    global H;
    global U_Long;

    
    %----------------------------------------------------------------------
    % input arguments
    %----------------------------------------------------------------------
    if isempty(varargin{1})
        error('lftfun is missing');
    else
        ext_lftfun = varargin{1};
    end;
    if nargin==1
        error('InputKnown is missing or bad input argument order');
    else
        ext_InputKnown = varargin{2};
    end;
    if nargin==2
        error('InitialConditions are missing or bad input argument order');
    else
        ext_InitialConditions = varargin{3};
    end;
    if nargin==3
        error('InputOptStart is missing or bad input argument order');
    else
        ext_InputOptStart = varargin{4};
    end;
    if nargin==4
        error('YToFollow is missing or bad input argument order');
    else
        ext_YToFollow = varargin{5};                                        % same struct as Input
    end;
    if nargin==5
        error('UObjective is missing or bad input argument order');
    else
        ext_UObjective = varargin{6};                                       % same struct as Input
    end;
    if (nargin==6)||isempty(varargin{7})
        ext_lftSolverOptions = lftSet();
    else
        ext_lftSolverOptions = varargin{7};
    end;
    if (nargin==7)||isempty(varargin{8})
        ext_lftOptimOptions = lftOptSet();
    else
        ext_lftOptimOptions = varargin{8};
    end;

 
    %% check and preparation
    num_u = size(ext_lftfun.LTI.B3,2);
    num_x = size(ext_lftfun.LTI.A,2);
    num_y = size(ext_lftfun.LTI.C3,1);

    num_u_o = size(ext_InputOptStart.Coefficients,1);                       % gestire num_u_o < num_u e ordine diverso da ISO e intrecci
    IndexInputToBeOptimized = [1,2];                                        % gestire num_u_o < num_u e ordine diverso da ISO e intrecci
    Input = ext_InputOptStart;                                              % gestire num_u_o < num_u e ordine diverso da ISO e intrecci
    
    MaxIter = ext_lftOptimOptions.MaxIter;
    EpsilonLambda = ext_lftOptimOptions.EpsilonLambda;
    MinNormGrad = ext_lftOptimOptions.MinNormGrad;
    RelTolX = ext_lftOptimOptions.RelTolX;
    NumberOfSteps = ext_lftOptimOptions.NumberOfSteps;
    ext_lftSolverOptions.NumberOfSteps = (NumberOfSteps + ext_lftSolverOptions.MaxOrder - 2) + 1;
    Input.Time = [0, ext_lftSolverOptions.MaxStep*ext_lftSolverOptions.NumberOfSteps];
    NumberOfData = (1+NumberOfSteps)*num_u_o;


    Q = eye(num_y*(1+NumberOfSteps));
    R = eye(num_u_o*(1+NumberOfSteps));
    
    A = zeros(num_x,num_x,1+NumberOfSteps);
    B = zeros(num_x,num_u_o,1+NumberOfSteps);
    dYdU = zeros(num_y*(1+NumberOfSteps),num_u_o*(1+NumberOfSteps));

    %----------------------------------------------------------------------
    % Optimization
    %----------------------------------------------------------------------
    disp('Continuous-Time Optimal Control task');
    disp('-----------------------------------------------------------------------------------------------------------------');
    disp('|  Iter  |     f(x)     |    LocalConditions    |         Method         |    NormOfStep    |   ControlWeight   |');
    disp('-----------------------------------------------------------------------------------------------------------------');
    formatstr = '    %3.0d   %13.6g   %20s %25s   %13.6g        %13.6g   '; 
    j_history = [];
    iter=1;
    delta_x = 0;
    while iter<=MaxIter
        CostFunction();
        j_history = cat(1,j_history,J);  
        [V,E] = eig(H);
        [lambda,LambdaIndex] = sort(diag(E));
        V = V(:,LambdaIndex);
        AbsLambda = abs(diag(E));
        delta_x_prev = delta_x;
        if lambda(1)>EpsilonLambda                                          % the problem is locally convex
            Method = 'Gauss-Newton';
            LocalConditions = 'Convex';
            [SortAbsLambda,~] = sort(AbsLambda);
            conditionNumber = abs(SortAbsLambda(end)/SortAbsLambda(1));
            h = 1;%/conditionNumber;
            delta_x = - (V*diag(1./lambda)*V')*G'*h;                        % move using standard Gauss-Newton algorithm
        else
            LambdaIndex = find(AbsLambda>EpsilonLambda);
            if length(LambdaIndex)==NumberOfData
                LocalConditions = 'Non-Convex';                             % the problem is locally non-convex and without singularity directions
                W = V;
            elseif ~isempty(LambdaIndex)
                LocalConditions = 'Singular';                               % the problem has singularity directions and convexity or concavity directions
                W = V(:,LambdaIndex);
            else
                LocalConditions = 'Singular (flat)';                        % the problem has nearly-zero curvature in any direction 
                W = eye(NumberOfData);
            end;
            gh = W'*G';                                                     % gh is the gradient projection on the non-null eigenspace
            if norm(G)>MinNormGrad                                          % the gradient has an important impact on the solution research
                [MaxGrad,iMaxGrad] = max(G);
                [MinGrad,iMinGrad] = min(G);
                [MaxGradNewBase,iMaxGradNewBase] = max(gh);
                [MinGradNewBase,iMinGradNewBase] = min(gh);
                [~,iMaximumGrad] = max(abs([MaxGrad,MinGrad,MaxGradNewBase,MinGradNewBase]));
                switch iMaximumGrad(1)
                    case 1
                        Method = 'Max Grad Dir';
                        GuessDirection = zeros(NumberOfData,1);
                        if G(iMaxGrad)>0
                            GuessDirection(iMaxGrad) = -1;
                        else
                            GuessDirection(iMaxGrad) = 1;
                        end;
                    case 2
                        Method = 'Max Grad Dir';
                        GuessDirection = zeros(NumberOfData,1);
                        if G(iMinGrad)>0
                            GuessDirection(iMinGrad) = -1;
                        else
                            GuessDirection(iMinGrad) = 1;
                        end;
                    case 3
                        Method = 'Max Grad Dir (eig. b.)';
                        if gh(iMaxGradNewBase)>0
                            GuessDirection = -W(:,iMaxGradNewBase);
                        else
                            GuessDirection = W(:,iMaxGradNewBase);
                        end;
                    case 4
                        Method = 'Max Grad Dir (eig. b.)';
                        if gh(iMinGradNewBase)>0
                            GuessDirection = -W(:,iMinGradNewBase);
                        else
                            GuessDirection = W(:,iMinGradNewBase);
                        end;
                end;
                h = 1;
                delta_x = GuessDirection*h;                                 % move following maximum absolute gradient direction and decreasing sense
                if norm(delta_x_prev+delta_x)<=norm(delta_x_prev)           % if the step from the actual point makes the next point closer to the previous one...
                    h = 0.5;
                    delta_x = GuessDirection*h;                             % ...move following maximum absolute gradient direction and decreasing sense with half versor magnitude                   
                end;
            else                                                            % the gradient nearly null in any direction
                Method = 'Max Conc Dir';            
                [MaxCoord,iMaxCoord] = max(V(1));
                [MinCoord,iMinCoord] = min(V(1));
                if abs(MaxCoord)>abs(MinCoord)
                    if G(iMaxCoord)<0
                        GuessDirection = V(:,1);
                    else
                        GuessDirection = -V(:,1);
                    end;
                else
                    if G(iMinCoord)>0
                        GuessDirection = V(:,1);
                    else
                        GuessDirection = -V(:,1);
                    end;
                end;
                h = 1;
                delta_x = GuessDirection*h;                                 % move following the maximum concavity direction using the sign of the gradient corresponding to the maximum absolute direction component
                if norm(delta_x_prev+delta_x)<=norm(delta_x_prev)           % if the step from the actual point makes the next point closer to the previous one...
                    if strcmp(LocalConditions, 'Singular (flat)')           % ...and the problem has nearly-zero curvature in any direction...
                        LocalConditions = 'Flat Region';                    % ...the problem is locally flat,... 
                        Method = 'Skip zone';                               % ...the research of the minimum stops to restart from an another point
                        normOfStep = 0;
                        h = 0;
                        iter_info = sprintf(formatstr,iter,J,LocalConditions,Method,normOfStep,h);
                        disp(iter_info);
                        break;
                    else                                                    % ...and the curvature has an important impact on the solution...
                        h = 0.5;
                        delta_x = GuessDirection*h;                         % ...move following the maximum concavity direction using the sign of the gradient corresponding to the maximum absolute direction component. The magnitude is half versor
                    end;
                end;
            end;
        end;
        normOfStep = norm(delta_x);
        x = U_Long;
        if (normOfStep/norm(x)<RelTolX)                                     % if the norm of the step to be used is too small with respect to the norm of the local solution...
            LocalConditions = 'Local Minimum Found';
            Method = 'Stop';
            normOfStep = 0;
            h = 0;
            iter_info = sprintf(formatstr,iter,J,LocalConditions,Method,normOfStep,h);
            disp(iter_info);
            break;                                                          % ...the research of the minimum stops because a local minimum has been found
        end;
        iter_info = sprintf(formatstr,iter,J,LocalConditions,Method,normOfStep,h);
        disp(iter_info);
        x = x + delta_x;                                                    % solution refresh
        
        
        
        % Unknown Input refresh
        for t=1:1+NumberOfSteps
            Input.Samples(t,:) = x((t-1)*num_u_o+1:t*num_u_o,1)';
        end;
        
        figure
        subplot 211
        plot(SensTime, Input.Samples(:,1), 'Marker', 'x', 'MarkerSize',10, 'MarkerEdgeColor','k');
        subplot 212
        plot(SensTime, Input.Samples(:,2),'Marker', 'x', 'MarkerSize',10, 'MarkerEdgeColor','k');        
        for uu=1:num_u_o
            Input.Coefficients(IndexInputToBeOptimized(uu),:) = polyfit(SensTime, Input.Samples(:,IndexInputToBeOptimized(uu)),ext_lftOptimOptions.PolynomialDegree);  
        end;
   
        iter = iter+1;
    end;
    disp('-----------------------------------------------------------------------------------------------------------------');
    Input.Time = SensTime;
    CostFun = J;

    
    

    
    %----------------------------------------------------------------------
    % CostFunction
    %----------------------------------------------------------------------
    function CostFunction()
        % LFT basic call
        ext_lftSolverOptions.Sensitivity = 0;
        [~, InternalSolution, CommonTerms] = lftSolver(ext_lftfun,Input,ext_InitialConditions,ext_lftSolverOptions);
        InternalSolution = [InternalSolution(1,:); InternalSolution(ext_lftSolverOptions.MaxOrder:end-1,:)];   
        CommonTerms.VCT1 = cat(1, CommonTerms.VCT1(1,:), CommonTerms.VCT1(ext_lftSolverOptions.MaxOrder:end-1,:));
        CommonTerms.VCT2 = cat(3, CommonTerms.VCT2(:,:,1), CommonTerms.VCT2(:,:,ext_lftSolverOptions.MaxOrder:end-1));
        CommonTerms.VCT3 = cat(1, CommonTerms.VCT3(1,:), CommonTerms.VCT3(ext_lftSolverOptions.MaxOrder:end-1,:));
        for cc=4:10
            eval(['CommonTerms.VCT',num2str(cc),' = cat(3, CommonTerms.VCT',num2str(cc),'(:,:,1), CommonTerms.VCT',num2str(cc),'(:,:,ext_lftSolverOptions.MaxOrder:end-1));']);
        end;
        % Input sensitivity routine
        SensOutput = [];
        ext_lftSolverOptions.CommonTerms = CommonTerms;
        SensTime = InternalSolution(:,1);
        for s=1:num_u_o
            ext_lftSolverOptions.Sensitivity = strcat('u', num2str(IndexInputToBeOptimized(s)));
            [~,is] = lftSolver(ext_lftfun,SensTime,ext_InitialConditions,ext_lftSolverOptions);
            SensOutput = cat(3,SensOutput,is);
        end;
        % A B matrix building (backwardEuler used for discretization)
        delta_t = 0.1*2/max(sort(abs(eig(CommonTerms.VCT10(1:num_x,1:num_x,1)))));
        A(:,:,1) = (eye(num_x)-delta_t*CommonTerms.VCT10(1:num_x,1:num_x,1))\eye(num_x);
        B(:,:,1) = A(:,:,1)*delta_t*mysqueeze(SensOutput(1,1+num_x+1:1+num_x+num_x,:))';
        for tt=2:1+NumberOfSteps
            delta_t = (SensTime(tt)-SensTime(tt-1));
            A(:,:,tt) = (eye(num_x)-delta_t*CommonTerms.VCT10(1:num_x,1:num_x,tt))\eye(num_x);
            B(:,:,tt) = A(:,:,tt)*delta_t*mysqueeze(SensOutput(tt,1+num_x+1:1+num_x+num_x,:))';
        end;
        % dYdU matrix building
        for j=1:1+NumberOfSteps
            dYdU((j-1)*num_y+1:j*num_y,(j-1)*num_u_o+1:j*num_u_o) = mysqueeze(SensOutput(j,end-num_y+1:end,:))';     
            T = B(:,:,j);
            for i=j+1:1+NumberOfSteps
                T = A(:,:,i)*T;
                dYdU((i-1)*num_y+1:i*num_y,(j-1)*num_u_o+1:j*num_u_o) = mysqueeze(CommonTerms.VCT10(num_x+1:end,1:num_x,i))*T;
            end;
        end;
        % U_SetPoint and Y_SetPoint evaluation
        U_SetPoint = zeros(1+NumberOfSteps, num_u_o);
        Y_SetPoint = zeros(1+NumberOfSteps, num_y);
        for jj=1:num_u_o
            U_SetPoint(:,jj) = polyval(ext_UObjective.Coefficients(jj,:),SensTime);
        end;
        for jj=1:num_y
            Y_SetPoint(:,jj) = polyval(ext_YToFollow.Coefficients(jj,:),SensTime);
        end;
        % CostFunction, Gradient, Hessian building
        U_T = InternalSolution(:,1+1:1+num_u_o)';
        U_Long = U_T(:);
        Y_T = InternalSolution(:,end-num_y+1:end)';
        U0_T = U_SetPoint';
        Y0_T = Y_SetPoint';
        G1 = (Y_T(:)-Y0_T(:))'*Q;                                           % [u1(0), u2(0)... u1(1), u2(1)...] is obtained by transposition and after (:)
        G2 = (U_Long-U0_T(:))'*R;
        % costruzione di G, H e J
        G=(1/NumberOfSteps)*(G2+G1*dYdU);                    
        H=(1/NumberOfSteps)*(R+dYdU'*Q*dYdU);
        J=(1/2)*(1/NumberOfSteps)*(G2*(U_Long-U0_T(:))+G1*(Y_T(:)-Y0_T(:)));
    end
    %----------------------------------------------------------------------





end








