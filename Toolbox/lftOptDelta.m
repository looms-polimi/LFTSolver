function [Delta_opt fval CostFunction gradient hessian ConditionNumber par_history, varargout] = lftOptDelta(varargin)



    %% input arguments
    if isempty(varargin{1})
        error('lftfun is missing');
    else
        ext_lftfun = varargin{1};
    end;
    num_u = size(ext_lftfun.LTI.B3,2);
    num_x = size(ext_lftfun.LTI.A,2);
    num_z = size(ext_lftfun.LTI.B1,2);
    num_y = size(ext_lftfun.LTI.C3,1);
    indToIdentify = find([ext_lftfun.DeltaSym{2:end,6}]);
    num_uncertain_par = length(indToIdentify);
        % rescaling conditions
        if ~isfield(ext_lftfun, 'ISO')
            Rescale = false;
        else
            if ~isfield(ext_lftfun.ISO, 'Nominal')
                Rescale = false;
            else
                if isrow(ext_lftfun.ISO.Nominal)
                    NominalVector = ext_lftfun.ISO.Nominal;
                else
                    NominalVector = ext_lftfun.ISO.Nominal';
                end;
                if isequal(NominalVector,ones(1,num_y))
                    Rescale = false;
                else
                    Rescale = true;
                    nonZero = find(NominalVector);
                    NominalMatrix = repmat(NominalVector(nonZero)',1,num_uncertain_par);
                end;
                
            end;
        end;
        % choise for Hessian and Gradient building
        buildingType = 0;
        if (num_y>1)&&(num_uncertain_par>1)&&(Rescale)
            buildingType = 1;
        end;
        if (num_y>1)&&(num_uncertain_par>1)&&(~Rescale)
            buildingType = 2;
        end;
        if (num_y>1)&&(num_uncertain_par==1)&&(Rescale)
            buildingType = 3;
        end;
        if (num_y>1)&&(num_uncertain_par==1)&&(~Rescale)
            buildingType = 4;
        end;
        if (num_y==1)&&(num_uncertain_par>1)
            buildingType = 5;
        end;
        if (num_y==1)&&(num_uncertain_par==1)
            buildingType = 6;
        end;
    if isempty(varargin{2})
        error('Input is missing');
    else
        ext_Input = varargin{2};
    end;
    if isempty(varargin{3})
        error('Initial conditions are missing');
    else
        ext_InitialConditions = varargin{3};
    end;
    SensitivityInitialConditions = struct('StateInitialConditions',zeros(num_x,1));
    if isempty(varargin{4})
        ext_LFTsolverOptions = lftSet('SolutionTimeSpan', ext_Input.Time);
    else
        ext_LFTsolverOptions = varargin{4};
        ext_LFTsolverOptions.SolutionTimeSpan = ext_Input.Time;
    end;
    if isempty(varargin{5})
        error('Sample are missing');
    else
        y_measured = varargin{5};
    end;
    if nargin == 5
        ext_LFToptimOptions = lftOptSet;
    else
        if isempty(varargin{6})
            ext_LFToptimOptions = lftOptSet;
        else
            ext_LFToptimOptions = varargin{6};
        end;
    end;
    
    
    %% optimization
    global J;
    global H;
    global g;
    par_history = [];
    ConditionNumber = [];
    gradient = [];
    hessian = [];
    CostFunction = [];
    options = optimset('Algorithm','trust-region-reflective',...
                       'OutputFcn', @myoutput,...
                       'TolFun', ext_LFToptimOptions.TolFun,...
                       'MaxIter',ext_LFToptimOptions.MaxIter,...
                       'Display',ext_LFToptimOptions.Display,...
                       'GradObj','on',...
                       'Hessian','on');
                   %                        'TolX',ext_LFToptimOptions.StepTolerance,...
    disp('LFT identification task');
    delta_vector_start = diag(ext_lftfun.DeltaVal);
    lb = [ext_lftfun.DeltaSym{[1+indToIdentify],1+6}]';
    ub = [ext_lftfun.DeltaSym{[1+indToIdentify],1+7}]';
    [delta_vector_opt, fval] = fmincon(@cost_function, delta_vector_start([ext_lftfun.DeltaSym{1+indToIdentify,2}]), [],[],[],[],lb,ub,[],options);
    extended_delta_vector_opt = delta_vector_start;
    for dd=1:num_uncertain_par
        extended_delta_vector_opt(ext_lftfun.DeltaSym{1+indToIdentify(dd),2}:ext_lftfun.DeltaSym{1+indToIdentify(dd),3}) = delta_vector_opt(dd);
    end;
    Delta_opt = diag(extended_delta_vector_opt);


    %% output function
    function stop = myoutput(delta_vector_opt,~,state)
        stop = false;
        if strcmp(state,'iter')
          par_history = cat(1, par_history, delta_vector_opt');
          ConditionNumber = cat(2, ConditionNumber, cond(H));
          gradient = cat(1, gradient, g);
          hessian = cat(3, hessian, H);
          CostFunction = cat(1, CostFunction, J);
        end
    end


    %% cost function
    function [J, g, H] = cost_function(Delta)
       global J;
       global H;
       global g;
       SensOutput = [];
       g = zeros(1,num_uncertain_par);
       J = 0;
       H = zeros(num_uncertain_par,num_uncertain_par);
       n_sample = 0; 
       extended_delta_vector = diag(ext_lftfun.DeltaVal);
       for d=1:num_uncertain_par
           extended_delta_vector(ext_lftfun.DeltaSym{1+indToIdentify(d),2}:ext_lftfun.DeltaSym{1+indToIdentify(d),3}) = Delta(d);
       end;
       ext_lftfun.DeltaVal = diag(extended_delta_vector);
       ext_LFTsolverOptions.Sensitivity = 0;
       [Output, InternalSolution, CommonTerms] = lftSolver(ext_lftfun,ext_Input,ext_InitialConditions,ext_LFTsolverOptions);
       SensInput = struct('Type', 'interpolated',...
                          'Time', InternalSolution(:,1),...
                          'Samples', InternalSolution(:,1+num_u+num_x+1:1+num_u+num_x+num_z));
       ext_LFTsolverOptions.CommonTerms = CommonTerms;
       for s = 1:num_uncertain_par
          ext_LFTsolverOptions.Sensitivity = strcat('p', num2str(indToIdentify(s)));
          SensOutput = cat(3,SensOutput,lftSolver(ext_lftfun,SensInput,SensitivityInitialConditions,ext_LFTsolverOptions));  
       end;
       for i=ext_LFToptimOptions.StartOptimSample:length(ext_Input.Time)  
           if isempty(strfind(num2str(y_measured(i,:)),'NaN'))
               n_sample = n_sample+1;
               % permute the vector SensOutput (containing the simulation 
               % time as second dimension) can jeopardize the efficency!
               switch buildingType
                   case 1
                       % num_y>1 AND num_uncertain_par>1 AND rescale
                       e = (y_measured(i,nonZero) - Output(i,end-num_y+nonZero))./NominalVector(nonZero);
                       if length(nonZero)>1
                           dydd = squeeze(SensOutput(i,end-num_y+nonZero,:))./NominalMatrix;
                       else
                           dydd = squeeze(SensOutput(i,end-num_y+nonZero,:))'./NominalMatrix;
                       end;
                   case 2
                       % num_y>1 AND num_uncertain_par>1 AND not_rescale
                       e = y_measured(i,:) - Output(i,end-num_y+1:end);
                       dydd = squeeze(SensOutput(i,end-num_y+1:end,:));
                   case 3
                       % num_y>1 AND num_uncertain_par=1 AND rescale
                       e = (y_measured(i,nonZero) - Output(i,end-num_y+nonZero))./NominalVector(nonZero);
                       dydd = (SensOutput(i,end-num_y+nonZero:end)./NominalVector(nonZero))';
                   case 4
                       % num_y>1 AND num_uncertain_par=1 AND not_rescale
                       e = y_measured(i,:) - Output(i,end-num_y+1:end);         %%%%%%%%%%%%%%% nonZero?????????????
                       dydd = SensOutput(i,end-num_y+1:end)';
                   case 5
                       % num_y=1 AND num_uncertain_par>1
                       e = y_measured(i,:) - Output(i,end);
                       dydd = (squeeze(SensOutput(i,end,:)))';
                   case 6
                       % num_y=1 AND num_uncertain_par=1
                       e = y_measured(i,:) - Output(i,end);
                       dydd = SensOutput(i,end);
                   otherwise
                       disp('internal LFToptimizer error: case not found for the building of the Hessian and the Gradient');
               end;
               H = H + dydd'*dydd;
               g = g - e*dydd;
               J = J + e*e';
          end;
       end;
       J = J*0.5/n_sample;
       g = g/n_sample;
       H = H/n_sample;
   end
end
