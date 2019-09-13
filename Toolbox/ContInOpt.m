function [y,x_history,j_history]  = ContInOpt(varargin)

    %----------------------------------------------------------------------
    % input argument
    %----------------------------------------------------------------------
    x = varargin{1};
    MaxIter = varargin{2};
    EpsilonLambda = varargin{3};
    RelTolX = varargin{4};
    MinNormGrad = 1e-6;
    NumberOfSamples = varargin{5};
    %----------------------------------------------------------------------

    
    
    disp('-----------------------------------------------------------------------------------------------------------------');
    disp('|  Iter  |     f(x)     |    LocalConditions    |         Method         |    NormOfStep    |   ControlWeight   |');
    disp('-----------------------------------------------------------------------------------------------------------------');
    formatstr = '    %3.0d   %13.6g   %20s %25s   %13.6g        %13.6g   '; 
    x_history = [];
    j_history = [];
    iter=1;
    delta_x = 0;
    while iter<=MaxIter

        J = CostFunction();
        G = Gradient();
        H = Hessian();

        x_history = cat(2,x_history,x);
        j_history = cat(1,j_history,J);  
        [V,D] = eig(H);
        [lambda,LambdaIndex] = sort(diag(D));
        V = V(:,LambdaIndex);
        AbsLambda = abs(diag(D));
        delta_x_prev = delta_x;
        if lambda(1)>EpsilonLambda                                          % the problem is locally convex
            Method = 'Gauss-Newton';
            LocalConditions = 'Convex';
            [SortAbsLambda,~] = sort(AbsLambda);
            conditionNumber = abs(SortAbsLambda(end)/SortAbsLambda(1));
            h = 1/conditionNumber;
            delta_x = - (V*diag(1./lambda)*V')*G'*h;                        % move using standard Gauss-Newton algorithm
        else
            LambdaIndex = find(AbsLambda>EpsilonLambda);
            if length(LambdaIndex)==NumberOfSamples
                LocalConditions = 'Non-Convex';                             % the problem is locally non-convex and without singularity directions
                W = V;
            elseif ~isempty(LambdaIndex)
                LocalConditions = 'Singular';                               % the problem has singularity directions and convexity or concavity directions
                W = V(:,LambdaIndex);
            else
                LocalConditions = 'Singular (flat)';                        % the problem has nearly-zero curvature in any direction 
                W = eye(NumberOfSamples);
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
                        GuessDirection = zeros(NumberOfSamples,1);
                        if G(iMaxGrad)>0
                            GuessDirection(iMaxGrad) = -1;
                        else
                            GuessDirection(iMaxGrad) = 1;
                        end;
                    case 2
                        Method = 'Max Grad Dir';
                        GuessDirection = zeros(NumberOfSamples,1);
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
        iter = iter+1;
    end;
    disp('-----------------------------------------------------------------------------------------------------------------');

    %%---------------------------------------------------------------------
    function y=CostFunction()
        y = 3*(1-x(1)).^2.*exp(-(x(1).^2) - (x(2)+1).^2)  - 10*(x(1)/5 - x(1).^3 - x(2).^5).*exp(-x(1).^2-x(2).^2) - 1/3*exp(-(x(1)+1).^2 - x(2).^2);
    end

    %%---------------------------------------------------------------------
    function y=Gradient()
        djdx1 =(exp(- (x(1) + 1)^2 - x(2)^2)*(2*x(1) + 2))/3 + 3*exp(- (x(2) + 1)^2 - x(1)^2)*(2*x(1) - 2) + exp(- x(1)^2 - x(2)^2)*(30*x(1)^2 - 2) - 6*x(1)*exp(- (x(2) + 1)^2 - x(1)^2)*(x(1) - 1)^2 - 2*x(1)*exp(- x(1)^2 - x(2)^2)*(10*x(1)^3 - 2*x(1) + 10*x(2)^5);
        djdx2 =(2*x(2)*exp(- (x(1) + 1)^2 - x(2)^2))/3 + 50*x(2)^4*exp(- x(1)^2 - x(2)^2) - 3*exp(- (x(2) + 1)^2 - x(1)^2)*(2*x(2) + 2)*(x(1) - 1)^2 - 2*x(2)*exp(- x(1)^2 - x(2)^2)*(10*x(1)^3 - 2*x(1) + 10*x(2)^5);
        y = [djdx1, djdx2];
    end

    %%---------------------------------------------------------------------
    function y=Hessian()
       djdx1x1 = (2*exp(- (x(1) + 1)^2 - x(2)^2))/3 + 6*exp(- (x(2) + 1)^2 - x(1)^2) - (exp(- (x(1) + 1)^2 - x(2)^2)*(2*x(1) + 2)^2)/3 - 2*exp(- x(1)^2 - x(2)^2)*(10*x(1)^3 - 2*x(1) + 10*x(2)^5) + 60*x(1)*exp(- x(1)^2 - x(2)^2) - 6*exp(- (x(2) + 1)^2 - x(1)^2)*(x(1) - 1)^2 - 12*x(1)*exp(- (x(2) + 1)^2 - x(1)^2)*(2*x(1) - 2) - 4*x(1)*exp(- x(1)^2 - x(2)^2)*(30*x(1)^2 - 2) + 4*x(1)^2*exp(- x(1)^2 - x(2)^2)*(10*x(1)^3 - 2*x(1) + 10*x(2)^5) + 12*x(1)^2*exp(- (x(2) + 1)^2 - x(1)^2)*(x(1) - 1)^2;
       djdx2x1 = 6*x(1)*exp(- (x(2) + 1)^2 - x(1)^2)*(2*x(2) + 2)*(x(1) - 1)^2 - 2*x(2)*exp(- x(1)^2 - x(2)^2)*(30*x(1)^2 - 2) - 100*x(1)*x(2)^4*exp(- x(1)^2 - x(2)^2) - 3*exp(- (x(2) + 1)^2 - x(1)^2)*(2*x(1) - 2)*(2*x(2) + 2) - (2*x(2)*exp(- (x(1) + 1)^2 - x(2)^2)*(2*x(1) + 2))/3 + 4*x(1)*x(2)*exp(- x(1)^2 - x(2)^2)*(10*x(1)^3 - 2*x(1) + 10*x(2)^5);
       djdx1x2 = djdx2x1;
       djdx2x2 = (2*exp(- (x(1) + 1)^2 - x(2)^2))/3 + 200*x(2)^3*exp(- x(1)^2 - x(2)^2) - 200*x(2)^5*exp(- x(1)^2 - x(2)^2) - (4*x(2)^2*exp(- (x(1) + 1)^2 - x(2)^2))/3 - 2*exp(- x(1)^2 - x(2)^2)*(10*x(1)^3 - 2*x(1) + 10*x(2)^5) - 6*exp(- (x(2) + 1)^2 - x(1)^2)*(x(1) - 1)^2 + 3*exp(- (x(2) + 1)^2 - x(1)^2)*(2*x(2) + 2)^2*(x(1) - 1)^2 + 4*x(2)^2*exp(- x(1)^2 - x(2)^2)*(10*x(1)^3 - 2*x(1) + 10*x(2)^5);
       y = [djdx1x1, djdx1x2
            djdx2x1, djdx2x2];
    end

y = x;



end




