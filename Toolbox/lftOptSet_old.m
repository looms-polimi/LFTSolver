function y = lftOptSet(varargin)
                LocalOptions = struct('TolFun', 1e-3,...
                                      'MaxIter',20,...
                                      'Display', 'off',...
                                      'StepTolerance', 1e-6,...
                                      'StartOptimSample', 1);
                for i=1:nargin/2
                    if ischar(varargin{i*2})
                        LocalValue = char(varargin(i*2));
                    else
                        LocalValue = varargin{i*2};
                    end;
                    LocalOptions.(char(varargin(i*2-1))) = LocalValue;
                end;
                y = LocalOptions;
end