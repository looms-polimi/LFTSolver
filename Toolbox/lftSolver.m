function [Solution, InternalSolution, CommonTerms, varargout] = lftSolver(varargin)



    %% input arguments
    %----------------------------------------------------------------------
    if isempty(varargin{1})
        error('lftfun is missing');
    else
        LTI = varargin{1}.LTI;
        DELTA = varargin{1}.DeltaVal;
    end;
    if nargin == 3
        LFTsolverOptions = lftSet;
    else
        LFTsolverOptions = varargin{4};
    end;
    Sensitivity = LFTsolverOptions.Sensitivity;
    if Sensitivity(1)=='u'
        if isempty(varargin{2})
            error('Time span is missing');
        else
            Time = varargin{2};
        end;
    else
        if isempty(varargin{2})
            error('Input is missing');
        else
            inputType = varargin{2}.Type;
            switch inputType
                case 'discontinuous'
                    if Sensitivity==0
                        Input_S = varargin{2}.Samples;
                        Time = varargin{2}.Time;
                    else
                        error('Wrong Type for Input argument for Sensitivity simulation');
                    end;
                case 'interpolated'
                    Input_S = varargin{2}.Samples;
                    Time = varargin{2}.Time;
                case 'continuous'
                    if Sensitivity==0
                        Input_C = varargin{2}.Coefficients;
                        MaxInputDegree = size(Input_C,2)-1;
                        Time = varargin{2}.Time;
                    else
                        error('Wrong Type for Input argument for Sensitivity simulation');
                    end;
                otherwise
                    error('Wrong Type for Input argument');
            end;
        end;
    end;
    if isempty(varargin{3})
        error('Initial conditions are missing');
    else
        switch Sensitivity(1)
            case 0
                StateIC = varargin{3}.StateInitialConditions;
            case 'u'
                StateIC = varargin{3}.SensInitialConditions(:,str2double(Sensitivity(2)));
            case 'p'
                StateIC = varargin{3}.StateInitialConditions;
        end;
    end;
    
    
    ShowIntTime = LFTsolverOptions.ShowIntTime;
    SensAlgorithm = LFTsolverOptions.SensAlgorithm;
    SolutionInterpMethod = LFTsolverOptions.SolutionInterpMethod;
    OversamplingMethod = LFTsolverOptions.OversamplingMethod;
    InitialStep = LFTsolverOptions.InitialStep;
    if isempty(LFTsolverOptions.MaxStep)
        MaxStep = Time(end);
    else
        MaxStep = LFTsolverOptions.MaxStep;
    end;
    NumberOfSteps = LFTsolverOptions.NumberOfSteps;
    CommonTerms = LFTsolverOptions.CommonTerms;
    

    
    %% cardinalities
    %---------------------------------------------------------------------- 
    num_u = size(LTI.B3,2);
    num_x = size(LTI.A,1);
    num_z = size(LTI.C1,1);
    num_o = size(LTI.C2,1);
    num_y = size(LTI.C3,1);
    if Sensitivity==0
        num_input = num_u;
    else
        switch Sensitivity(1)
            case 'p'
                delta_index = str2double(Sensitivity(2:end));
                num_input = num_z;
                ind_start_z = cell2mat(varargin{1}.DeltaSym(1+delta_index,2));
                ind_stop_z = cell2mat(varargin{1}.DeltaSym(1+delta_index,3));
            case 'u'
                u_index = str2double(Sensitivity(2));
            otherwise
                error('Wrong value for "Sensitivity"');
        end;
    end;
    

    
    %% global variables
    %----------------------------------------------------------------------
    global u_step;
    global theta_step;
    global dertheta_step;
    global ind_step;
    global VCT1;
    global VCT2;
    global VCT3;
    global VCT4;
    global VCT5;
    global VCT6;
    global VCT7;
    global VCT4_step;
    global VCT5_step;
    global VCT6_step;
    global VCT7_step;
    global I_step;
    global timeInterpStartIndex;
    global stop;

    
    
    %% simulation
    %----------------------------------------------------------------------
    stop = 0;
    ind_step = 1;
    if Sensitivity==0
        % CommonTerms initialization
        FCT1 = LTI.B1*DELTA;
        FCT2 = (LTI.D11*DELTA - eye(num_z));
        FCT3 = LTI.D21*DELTA;
        FCT4 = LTI.D31*DELTA;
        switch inputType
            case 'discontinuous'
                VCT3 = Input_S(1,:);
            case 'interpolated'
                VCT3 = Input_S(1,:);
                data1 = Input_S;
                data1_sparsity = [];
                data2 = [];
            otherwise % continuous
                VCT3 = Input_C(:,end)';
        end;
        if isfield(varargin{3}, 'AlgebraicStartingValues')&&(~isempty(varargin{3}.AlgebraicStartingValues))
            AlgebraicStartingValues = varargin{3}.AlgebraicStartingValues;
        else
            AlgebraicStartingValues = [(eye(num_z)-LTI.D11*DELTA)\(LTI.C1*StateIC+LTI.D13*VCT3');
                                           LTI.C2*StateIC+LTI.D23*VCT3'];
        end;
        AlgebraicIC = SearchConsistentAlgebraicIC();
        VCT1 = varargin{1}.Theta(AlgebraicIC(num_z+1:num_z+num_o))';        
        VCT2 = varargin{1}.dThetadOmega(AlgebraicIC(num_z+1:num_z+num_o));
        VCT4 = LTI.B2*VCT2;
        VCT5 = LTI.D12*VCT2;
        VCT6 = (LTI.D22*VCT2-eye(num_o));
        VCT7 = LTI.D32*VCT2;
        % Initial conditions
        x = [StateIC;AlgebraicIC]';
        YIC = ( LTI.C3*StateIC + FCT4*AlgebraicIC(1:num_z) + LTI.D32*varargin{1}.Theta(AlgebraicIC(num_z+1:num_z+num_o)) + LTI.D33*VCT3');
        InternalSolution = [Time(1);VCT3';StateIC;AlgebraicIC; varargin{1}.Theta(AlgebraicIC(num_z+1:num_z+num_o));YIC]';
        % DAE solver setup
        M = zeros(num_x+num_z+num_o);
        M(1:num_x,1:num_x) = eye(num_x);
        options = odeset('Mass', M,...
                         'MassSingular', 'yes',...
                         'Jacobian', @(t,x)fjac_DAE(t,x),...
                         'RelTol', LFTsolverOptions.RelTol,...
                         'AbsTol', LFTsolverOptions.AbsTol,...
                         'MaxOrder', LFTsolverOptions.MaxOrder,...
                         'BDF', LFTsolverOptions.BDF,...
                         'InitialStep', InitialStep,...
                         'MaxStep', MaxStep,...
                         'OutputFcn', @OdeOutput,...
                         'Events', @StopSim);
        % Simulation routine
        index = 1;
        while index<length(Time)&&(~stop)
            X0 = x(end,:)';
            ind_start = index;
            switch inputType
                case 'discontinuous'
                    U0 = Input_S(ind_start,:)';
                    while isequal(Input_S(index,:)',U0)&&(index<length(Time))
                        index = index + 1;
                    end
                case 'interpolated'
                    U0 = Input_S;
                    index = length(Time);
                otherwise % continuous
                    U0 = Input_C;
                    index = length(Time);
            end;
            ind_stop = index;
            Tspan=[Time(ind_start),Time(ind_stop)];
            [t_ode15s,x] = ode15s(@(t,x)LFTmodel(t,x),Tspan,X0,options);
            y = zeros(length(t_ode15s),num_y);
            for i=2:length(t_ode15s)
                    y(i,:) = ( LTI.C3*x(i,1:num_x)' + FCT4*x(i,num_x+1:num_x+num_z)' + LTI.D32*VCT1(ind_step-length(t_ode15s)+i,:)' + LTI.D33*VCT3(ind_step-length(t_ode15s)+i,:)' )';
            end;  
            InternalSolution = cat(1, InternalSolution, [t_ode15s(2:end), VCT3(ind_step-length(t_ode15s)+2:ind_step,:), x(2:end,:), VCT1(ind_step-length(t_ode15s)+2:ind_step,:), y(2:end,:)]);
        end;
%% obsolete
%         % Variable common terms (without mem preallocation)
%         VCT8 = [];
%         VCT9_B = [];
%         VCT10 = [];
%         for i = 1:ind_step
%             VCT8 = cat(3, VCT8, ([-FCT2, -VCT5(:,:,i)
%                                   -FCT3, -VCT6(:,:,i)]\[LTI.C1, LTI.D13, LTI.D11
%                                                         LTI.C2, LTI.D23, LTI.D21]));
%             VCT9_B = cat(3, VCT9_B, ([LTI.B1
%                                       LTI.D31]+[FCT1, VCT4(:,:,i)
%                                                 FCT4, VCT7(:,:,i)]*VCT8(:,num_x+num_input+1:end,i)));
%             VCT10 = cat(3, VCT10, ([LTI.A,  LTI.B3
%                                     LTI.C3, LTI.D33]+[FCT1, VCT4(:,:,i)
%                                                       FCT4, VCT7(:,:,i)]*VCT8(:,1:num_x+num_input,i)));    
%         end;
%%
        % Variable common terms (with mem preallocation)
        VCT8 = zeros(num_z+num_o,num_x+num_u+num_z,ind_step);
        VCT9_B = zeros(num_x+num_y,num_z,ind_step);
        VCT10 = zeros(num_x+num_y,num_x+num_u,ind_step);
        for i = 1:ind_step
            VCT8(:,:,i) = [-FCT2, -VCT5(:,:,i)
                           -FCT3, -VCT6(:,:,i)]\[LTI.C1, LTI.D13, LTI.D11
                                                 LTI.C2, LTI.D23, LTI.D21];
            VCT9_B(:,:,i) = [LTI.B1
                             LTI.D31]+[FCT1, VCT4(:,:,i)
                                       FCT4, VCT7(:,:,i)]*VCT8(:,num_x+num_input+1:end,i);
            VCT10(:,:,i) = [LTI.A,  LTI.B3
                            LTI.C3, LTI.D33]+[FCT1, VCT4(:,:,i)
                                              FCT4, VCT7(:,:,i)]*VCT8(:,1:num_x+num_input,i);    
        end
        VCT9 = cat(2, VCT10(:,1:num_x,:), VCT9_B);
        CommonTerms = struct('FCT1', FCT1,...
                             'FCT2', FCT2,...
                             'FCT3', FCT3,...
                             'FCT4', FCT4,...
                             'VCT1', VCT1,...
                             'VCT2', VCT2,...
                             'VCT3', VCT3,...
                             'VCT4', VCT4,...
                             'VCT5', VCT5,...
                             'VCT6', VCT6,...
                             'VCT7', VCT7,...
                             'VCT8', VCT8,...
                             'VCT9', VCT9,...
                             'VCT10', VCT10);   
    else
        switch Sensitivity(1)
            case 'p'
                data1_sparsity = varargin{1}.Sparsity.VCT9{1,1};
                data2_sparsity = varargin{1}.Sparsity.VCT9{1,1+delta_index};
                data4_sparsity = varargin{1}.Sparsity.VCT9{2,1+delta_index};
                data5_sparsity = [];
                data1 = CommonTerms.VCT9(1:num_x,1:num_x,:);
                if data2_sparsity.NumNonZeroElements==0
                    data2 = [];
                else
                    data2 = CommonTerms.VCT9(1:num_x,num_x+ind_start_z:num_x+ind_stop_z,:);
                end;
                data3 = CommonTerms.VCT9(num_x+1:end,1:num_x,:);
                data4 = CommonTerms.VCT9(num_x+1:end,num_x+ind_start_z:num_x+ind_stop_z,:);
                data5 = Input_S(:,ind_start_z:ind_stop_z);
                InternalSolution = [Time, zeros(length(Time),num_x+num_y)];
                switch SensAlgorithm
                    case 'backwardEuler'
                        InternalSolution(1,2:1+num_x) = StateIC';
                        if data2_sparsity.NumNonZeroElements==0
                            if data4_sparsity.NumNonZeroElements==0
                                InternalSolution(1,1+num_x+1:end) = data3(:,:,1)*StateIC;
                                for i=2:length(Time)
                                    h = Time(i)-Time(i-1);
                                    InternalSolution(i,2:1+num_x) = (eye(num_x)-h*data1(:,:,i))\(InternalSolution(i-1,2:1+num_x)');                                
                                    InternalSolution(i,1+num_x+1:end) = data3(:,:,i)*InternalSolution(i,2:1+num_x)';
                                end;
                            else
                                InternalSolution(1,1+num_x+1:end) = data3(:,:,1)*StateIC +...
                                                                    data4(:,:,1)*data5(1,:)';
                                for i=2:length(Time)
                                    h = Time(i)-Time(i-1);
                                    InternalSolution(i,2:1+num_x) = (eye(num_x)-h*data1(:,:,i))\(InternalSolution(i-1,2:1+num_x)');                                
                                    InternalSolution(i,1+num_x+1:end) = data3(:,:,i)*InternalSolution(i,2:1+num_x)' +...
                                                                        data4(:,:,i)*data5(i,:)';
                                end;
                            end;
                        else
                            if data4_sparsity.NumNonZeroElements==0
                                InternalSolution(1,1+num_x+1:end) = data3(:,:,1)*StateIC;
                                for i=2:length(Time)
                                    h = Time(i)-Time(i-1);
                                    InternalSolution(i,2:1+num_x) = (eye(num_x)-h*data1(:,:,i))\(InternalSolution(i-1,2:1+num_x)' +...
                                                                    h*data2(:,:,i)*data5(i,:)');
                                    InternalSolution(i,1+num_x+1:end) = data3(:,:,i)*InternalSolution(i,2:1+num_x)';
                                end;
                            else
                                InternalSolution(1,1+num_x+1:end) = data3(:,:,1)*StateIC +...
                                                                    data4(:,:,1)*data5(1,:)';
                                for i=2:length(Time)
                                    h = Time(i)-Time(i-1);
                                    InternalSolution(i,2:1+num_x) = (eye(num_x)-h*data1(:,:,i))\(InternalSolution(i-1,2:1+num_x)' +...
                                                                    h*data2(:,:,i)*data5(i,:)');
                                    InternalSolution(i,1+num_x+1:end) = data3(:,:,i)*InternalSolution(i,2:1+num_x)' +...
                                                                        data4(:,:,i)*data5(i,:)';
                                end;
                            end;
                        end;
                    case 'ode15s'
                        options = odeset('Jacobian', @(t,x)fjac_DAE(t,x),...
                                         'RelTol', LFTsolverOptions.RelTol,...
                                         'AbsTol', LFTsolverOptions.AbsTol,...
                                         'InitialStep', InitialStep);     
                        [t_ode15s,x] = ode15s(@(t,x)LFTmodel(t,x),[Time(1),Time(end)],StateIC,options);    
                        for i=1:num_x
                            InternalSolution(:,1+i) = interp1(t_ode15s,x(:,i),Time,OversamplingMethod);
                        end;
                        if data4_sparsity.NumNonZeroElements==0
                            for i=1:length(Time)
                                InternalSolution(i,1+num_x+1:end) = data3(:,:,i)*InternalSolution(i,2:1+num_x)';
                            end;
                        else
                            for i=1:length(Time)
                                InternalSolution(i,1+num_x+1:end) = data3(:,:,i)*InternalSolution(i,2:1+num_x)' +...
                                                                    data4(:,:,i)*data5(i,:)';
                            end;
                        end;
                    otherwise
                        error('SensAlgorithm must be "backwardEuler" or "ode15s"');
                end;
            case 'u'
                data1_sparsity = varargin{1}.Sparsity.VCT10{1,1};
                data2_sparsity = varargin{1}.Sparsity.VCT10{1,1+u_index};
                data4_sparsity = varargin{1}.Sparsity.VCT10{2,1+u_index};
                data1 = CommonTerms.VCT10(1:num_x,1:num_x,:);
                if data2_sparsity.NumNonZeroElements==0
                    data2 = [];
                else
                    data2 = CommonTerms.VCT10(1:num_x,num_x+u_index,:);
                end;
                data3 = CommonTerms.VCT10(num_x+1:end,1:num_x,:);
                data4 = CommonTerms.VCT10(num_x+1:end,num_x+u_index,:);
                data5 = [];
                InternalSolution = [Time, zeros(length(Time),num_x+num_x+num_y)];
                switch SensAlgorithm
                    case 'backwardEuler'
                        InternalSolution(1,2:1+num_x) = StateIC';
                        if data2_sparsity.NumNonZeroElements==0
                            if data4_sparsity.NumNonZeroElements==0
                                InternalSolution(1,1+num_x+1:1+num_x+num_x) = data1(:,:,1)*StateIC;
                                InternalSolution(1,1+num_x+num_x+1:end) = data3(:,:,1)*StateIC;
                                for i=2:length(Time)
                                    h = Time(i)-Time(i-1);
                                    hg11 = h*data1(:,:,i);
                                    InternalSolution(i,2:1+num_x) = (eye(num_x)-hg11)\(InternalSolution(i-1,2:1+num_x)');
                                    InternalSolution(i,1+num_x+1:1+num_x+num_x) = data1(:,:,i)*InternalSolution(i,2:1+num_x)';
                                    InternalSolution(i,1+num_x+num_x+1:end) = data3(:,:,i)*InternalSolution(i,2:1+num_x)';                                          
                                end;                                
                            else
                                InternalSolution(1,1+num_x+1:1+num_x+num_x) = data1(:,:,1)*StateIC;
                                InternalSolution(1,1+num_x+num_x+1:end) = data3(:,:,1)*StateIC +...
                                                                          data4(:,:,1);
                                for i=2:length(Time)
                                    h = Time(i)-Time(i-1);
                                    hg11 = h*data1(:,:,i);
                                    InternalSolution(i,2:1+num_x) = (eye(num_x)-hg11)\(InternalSolution(i-1,2:1+num_x)');
                                    InternalSolution(i,1+num_x+1:1+num_x+num_x) = data1(:,:,i)*InternalSolution(i,2:1+num_x)';
                                    InternalSolution(i,1+num_x+num_x+1:end) = data3(:,:,i)*InternalSolution(i,2:1+num_x)' +...
                                                                              data4(:,:,i);                                          
                                end;
                            end;
                        else
                            if data4_sparsity.NumNonZeroElements==0
                                InternalSolution(1,1+num_x+1:1+num_x+num_x) = data1(:,:,1)*StateIC +...
                                                                              data2(:,:,1);
                                InternalSolution(1,1+num_x+num_x+1:end) = data3(:,:,1)*StateIC;
                                for i=2:length(Time)
                                    h = Time(i)-Time(i-1);
                                    hg11 = h*data1(:,:,i);
                                    hg12 = h*data2(:,:,i);
                                    InternalSolution(i,2:1+num_x) = (eye(num_x)-hg11)\(InternalSolution(i-1,2:1+num_x)' + hg12);
                                    InternalSolution(i,1+num_x+1:1+num_x+num_x) = data1(:,:,i)*InternalSolution(i,2:1+num_x)' +...
                                                                                  data2(:,:,i);
                                    InternalSolution(i,1+num_x+num_x+1:end) = data3(:,:,i)*InternalSolution(i,2:1+num_x)';                                          
                                end;
                            else
                                InternalSolution(1,1+num_x+1:1+num_x+num_x) = data1(:,:,1)*StateIC +...
                                                                              data2(:,:,1);
                                InternalSolution(1,1+num_x+num_x+1:end) = data3(:,:,1)*StateIC +...
                                                                          data4(:,:,1);
                                for i=2:length(Time)
                                    h = Time(i)-Time(i-1);
                                    hg11 = h*data1(:,:,i);
                                    hg12 = h*data2(:,:,i);
                                    InternalSolution(i,2:1+num_x) = (eye(num_x)-hg11)\(InternalSolution(i-1,2:1+num_x)' + hg12);
                                    InternalSolution(i,1+num_x+1:1+num_x+num_x) = data1(:,:,i)*InternalSolution(i,2:1+num_x)' +...
                                                                                  data2(:,:,i);
                                    InternalSolution(i,1+num_x+num_x+1:end) = data3(:,:,i)*InternalSolution(i,2:1+num_x)' +...
                                                                              data4(:,:,i);                                          
                                end;
                            end;
                        end;
                    case 'ode15s'
                        options = odeset('Jacobian', @(t,x)fjac_DAE(t,x),...
                                         'RelTol', LFTsolverOptions.RelTol,...
                                         'AbsTol', LFTsolverOptions.AbsTol,...
                                         'InitialStep', InitialStep);     
                        [t_ode15s,x] = ode15s(@(t,x)LFTmodel(t,x),[Time(1),Time(end)],StateIC,options);    
                        for i=1:num_x
                            InternalSolution(:,1+i) = interp1(t_ode15s,x(:,i),Time,OversamplingMethod);
                        end;
                        if data2_sparsity.NumNonZeroElements==0
                            if data4_sparsity.NumNonZeroElements==0
                                for i=1:length(Time)
                                    InternalSolution(i,1+num_x+1:1+num_x+num_x) = data1(:,:,i)*InternalSolution(i,2:1+num_x)';
                                    InternalSolution(i,1+num_x+num_x+1:end) = data3(:,:,i)*InternalSolution(i,2:1+num_x)';
                                end;
                            else
                                for i=1:length(Time)
                                    InternalSolution(i,1+num_x+1:1+num_x+num_x) = data1(:,:,i)*InternalSolution(i,2:1+num_x)';
                                    InternalSolution(i,1+num_x+num_x+1:end) = data3(:,:,i)*InternalSolution(i,2:1+num_x)' +...
                                                                              data4(:,:,i);
                                end;
                            end;
                        else
                            if data4_sparsity.NumNonZeroElements==0
                                for i=1:length(Time)
                                    InternalSolution(i,1+num_x+1:1+num_x+num_x) = data1(:,:,i)*InternalSolution(i,2:1+num_x)' +...
                                                                                  data2(:,:,i);
                                    InternalSolution(i,1+num_x+num_x+1:end) = data3(:,:,i)*InternalSolution(i,2:1+num_x)';
                                end;
                            else
                                for i=1:length(Time)
                                    InternalSolution(i,1+num_x+1:1+num_x+num_x) = data1(:,:,i)*InternalSolution(i,2:1+num_x)' +...
                                                                                  data2(:,:,i);
                                    InternalSolution(i,1+num_x+num_x+1:end) = data3(:,:,i)*InternalSolution(i,2:1+num_x)' +...
                                                                              data4(:,:,i);
                                end;
                            end;
                        end;
                    otherwise
                        error('SensAlgorithm must be "backwardEuler" or "ode15s"');
                end;
        end;
    end;
    if isempty(LFTsolverOptions.SolutionTimeSpan)
        Solution = [];
    else
        if ~isempty(NumberOfSteps)
            warning('SolutionTimeSpan ignored because of earlier termination required');
            Solution = [];
        else
            SolutionTimeSpan = LFTsolverOptions.SolutionTimeSpan;
            Solution = zeros(length(SolutionTimeSpan),1+num_y);
            Solution(:,1) = SolutionTimeSpan;
            for i=1:num_y
                Solution(:,1+i) = interp1(InternalSolution(:,1),InternalSolution(:,end-num_y+i),SolutionTimeSpan',SolutionInterpMethod);
            end;
        end;
    end;
    
    

            %% SearchConsistentAlgebraicInitialConditions
            %-----------------------------------------------------------------------------------
            function y = SearchConsistentAlgebraicIC()
                FSolveOptions = optimset('Display','off',...
                                         'Jacobian','on',...
                                         'TolFun',1e-12,...
                                         'TolX',1e-12);
                y = fsolve(@g0,AlgebraicStartingValues,FSolveOptions);
            end
     
        
        
            %% g0
            %-----------------------------------------------------------------------------------
            function [y,J]=g0(AlgebraicVariables)
                y(1:num_z)=LTI.C1*StateIC+FCT2*AlgebraicVariables(1:num_z)+LTI.D12*varargin{1}.Theta(AlgebraicVariables(num_z+1:num_z+num_o))+LTI.D13*VCT3';
                y(num_z+1:num_z+num_o)=LTI.C2*StateIC+FCT3*AlgebraicVariables(1:num_z)+LTI.D22*varargin{1}.Theta(AlgebraicVariables(num_z+1:num_z+num_o))+LTI.D23*VCT3'-AlgebraicVariables(num_z+1:num_z+num_o);
                J = [FCT2,    LTI.D12*varargin{1}.dThetadOmega(AlgebraicVariables(num_z+1:num_z+num_o))
                     FCT3,    LTI.D22*varargin{1}.dThetadOmega(AlgebraicVariables(num_z+1:num_z+num_o))-eye(num_o)];
            end 

        
            
            %% OdeOutput
            %-----------------------------------------------------------------------------------
            function status = OdeOutput(~,~,flag)
                if ~(strcmp(flag,'init')||strcmp(flag,'done'))
                    VCT1 = cat(1, VCT1, theta_step');
                    VCT2 = cat(3, VCT2, dertheta_step);
                    VCT3 = cat(1, VCT3, u_step');
                    VCT4 = cat(3, VCT4, VCT4_step);
                    VCT5 = cat(3, VCT5, VCT5_step);
                    VCT6 = cat(3, VCT6, VCT6_step);
                    VCT7 = cat(3, VCT7, VCT7_step);
                    ind_step = ind_step + 1;
                    status = 0;
                end;
            end
  
        
          
            %% LFTmodel
            %-----------------------------------------------------------------------------------
            function f = LFTmodel(t,x)
                if Sensitivity==0
                    x_local = x(1:num_x);
                    switch inputType
                        case 'discontinuous'
                            u_step = U0;
                        case 'interpolated'
                            u_local = timeInterp(t);
                            u_step = u_local{1}';
                        otherwise % continuous
                            u_local = U0(:,end);
                            tj = 1;
                            for jj=1:MaxInputDegree
                                tj = tj*t;
                                u_local = u_local + U0(:,end-jj)*tj;
                            end;
                            u_step = u_local;
                    end;
                    z_local = x(num_x+1:num_x+num_z);
                    omega_local = x(num_x+num_z+1:num_x+num_z+num_o);
                    theta_step = varargin{1}.Theta(omega_local);
                    dertheta_step = varargin{1}.dThetadOmega(omega_local);
                    VCT4_step = LTI.B2*dertheta_step;
                    VCT5_step = LTI.D12*dertheta_step;
                    VCT6_step = (LTI.D22*dertheta_step-eye(num_o));
                    VCT7_step = LTI.D32*dertheta_step;
                    f = [LTI.A*x_local + FCT1*z_local + LTI.B2*theta_step + LTI.B3*u_step
                         LTI.C1*x_local + FCT2*z_local + LTI.D12*theta_step + LTI.D13*u_step
                         LTI.C2*x_local + FCT3*z_local + LTI.D22*theta_step - eye(num_o)*omega_local + LTI.D23*u_step];                     
                elseif Sensitivity(1)=='p'
                    I_step = timeInterp(t);
                    f = I_step{1}*x + I_step{2}*I_step{3}';
                else
                    I_step = timeInterp(t);
                    f = I_step{1}*x + I_step{2};
                end;
                if ShowIntTime
                    t
                end;
            end
        
            
                     
            %% fjac_DAE
            %-----------------------------------------------------------------------------------        
            function y = fjac_DAE(~,~)
                if Sensitivity==0
                    y = [LTI.A,   FCT1,  VCT4_step
                         LTI.C1,  FCT2,  VCT5_step
                         LTI.C2,  FCT3,  VCT6_step];
                else
                    y = I_step{1};
                end;
            end
        
        
        
            %% timeInterp
            %-----------------------------------------------------------------------------------
            function y = timeInterp(t)
                if t==0
                    timeInterpStartIndex = 1;
                end;
                half_index = floor(length(Time)/2);
                half_time = Time(half_index);
                if t>=Time(timeInterpStartIndex)
                    if t>=half_time
                        [a,b] = timeFind(Time(timeInterpStartIndex:end),t,0);
                    else
                        [a,b] = timeFind(Time(timeInterpStartIndex:half_index),t,0);
                    end;
                    a = timeInterpStartIndex - 1 + a;
                    b = timeInterpStartIndex - 1 + b;
                else
                    if t>=half_time
                        [a,b] = timeFind(Time(half_index:timeInterpStartIndex),t,0);
                        a = half_index - 1 + a;
                        b = half_index - 1 + b;
                    else
                        [a,b] = timeFind(Time(1:timeInterpStartIndex),t,0);
                    end;
                end;
                q = (t-Time(a))/(Time(b)-Time(a));
                if size(data1,3)==1
                    y1 = zeros(size(data1(1,:)));
                    for ii=1:size(data1,2)
                        y1(ii)=data1(a,ii) + (data1(b,ii)-data1(a,ii))*q;
                    end;
                else
                    y1 = zeros(size(data1(:,:,1)));
                    if isempty(data1_sparsity)
                        for jj=1:size(data1,2)
                            for ii=1:size(data1,1)
                                y1(ii,jj)=data1(ii,jj,a) + (data1(ii,jj,b)-data1(ii,jj,a))*q;
                            end;
                        end;
                    else
                        for ii=1:data1_sparsity.NumNonZeroElements
                            y1(data1_sparsity.row(ii),data1_sparsity.col(ii))=data1(data1_sparsity.row(ii),data1_sparsity.col(ii),a)...
                            + (data1(data1_sparsity.row(ii),data1_sparsity.col(ii),b)-data1(data1_sparsity.row(ii),data1_sparsity.col(ii),a))*q;
                        end;
                    end;
                end;
                if ~isempty(data2)
                    if size(data2,3)==1
                        y2 = zeros(size(data2(1,:)));
                        for ii=1:size(data2,2)
                            y2(ii)=data2(a,ii) + (data2(b,ii)-data2(a,ii))*q;
                        end;
                    else
                        y2 = zeros(size(data2(:,:,1)));
                        if isempty(data2_sparsity)
                            for jj=1:size(data2,2)
                                for ii=1:size(data2,1)
                                    y2(ii,jj)=data2(ii,jj,a) + (data2(ii,jj,b)-data2(ii,jj,a))*q;
                                end;
                            end;
                        else
                            for ii=1:data2_sparsity.NumNonZeroElements
                                y2(data2_sparsity.row(ii),data2_sparsity.col(ii))=data2(data2_sparsity.row(ii),data2_sparsity.col(ii),a)...
                                + (data2(data2_sparsity.row(ii),data2_sparsity.col(ii),b)-data2(data2_sparsity.row(ii),data2_sparsity.col(ii),a))*q;
                            end;
                        end;
                    end;
                    if ~isempty(data5)
                        if size(data5,3)==1
                            y3 = zeros(size(data5(1,:)));
                            for ii=1:size(data5,2)
                                y3(ii)=data5(a,ii) + (data5(b,ii)-data5(a,ii))*q;
                            end;
                        else
                            y3 = zeros(size(data5(:,:,1)));
                            if isempty(data5_sparsity)
                                for jj=1:size(data5,2)
                                    for ii=1:size(data5,1)
                                        y3(ii,jj)=data5(ii,jj,a) + (data5(ii,jj,b)-data5(ii,jj,a))*q;
                                    end;
                                end;
                            else
                                for ii=1:data5_sparsity.NumNonZeroElements
                                    y3(data5_sparsity.row(ii),data5_sparsity.col(ii))=data5(data5_sparsity.row(ii),data5_sparsity.col(ii),a)...
                                    + (data5(data5_sparsity.row(ii),data5_sparsity.col(ii),b)-data5(data5_sparsity.row(ii),data5_sparsity.col(ii),a))*q;
                                end;
                            end;
                        end;
                        y = {y1,y2,y3};
                    else
                        y = {y1,y2,0};
                    end;
                else
                    y = {y1,0,0};
                end;
                timeInterpStartIndex = a;
            end
        
        
        
            %% StopSim
            %-----------------------------------------------------------------------------------
            function [ssval,isterminal,ssdir] = StopSim(~,~)
                if isempty(NumberOfSteps)
                    isterminal = 0;
                else 
                    isterminal = 1;
                end; 
                ssval = NumberOfSteps - ind_step;
                ssdir = 0;
                if ind_step == NumberOfSteps
                    stop = 1;
                end;
            end

        
    clear global u_step;
    clear global theta_step;
    clear global dertheta_step;
    clear global ind_step;
    clear global VCT1;
    clear global VCT2;
    clear global VCT3;
    clear global VCT4;
    clear global VCT5;
    clear global VCT6;
    clear global VCT7;
    clear global VCT4_step;
    clear global VCT5_step;
    clear global VCT6_step;
    clear global VCT7_step;
    clear global I_step;
    clear global timeInterpStartIndex;

    
    
end


