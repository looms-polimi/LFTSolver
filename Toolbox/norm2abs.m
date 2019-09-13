function [par_ValAbs_0, par_ValAbs_opt] = norm2abs(ext_lftfun, ext_DELTA_0, ext_DELTA_opt)

    num_uncertain_par = size(ext_lftfun.DeltaSym,1)-1;
    DELTA_0 = zeros(num_uncertain_par,1);
    DELTA_opt = zeros(num_uncertain_par,1);
    if isequal([ext_lftfun.DeltaSym{2:end,2}],[ext_lftfun.DeltaSym{2:end,3}])
        DELTA_0 = diag(ext_DELTA_0);
        DELTA_opt = diag(ext_DELTA_opt);
    else
        for dd=1:num_uncertain_par
            DELTA_0(dd,1) = ext_DELTA_0(ext_lftfun.DeltaSym{1+dd,2},ext_lftfun.DeltaSym{1+dd,2});
            DELTA_opt(dd,1) = ext_DELTA_opt(ext_lftfun.DeltaSym{1+dd,2},ext_lftfun.DeltaSym{1+dd,2});
        end;
    end;
    delta_min = [ext_lftfun.DeltaSym{2:end,4}];
    delta_max = [ext_lftfun.DeltaSym{2:end,5}];
    par_min = zeros(num_uncertain_par,1);
    D = zeros(num_uncertain_par,1);
    for ii=1:num_uncertain_par
        par_min(ii) = delta_min(ii);
        D(ii) = (delta_max(ii)-delta_min(ii));
    end;
    par_ValAbs_0 = par_min+0.5*diag(DELTA_0+ones(num_uncertain_par,1))*D;
    par_ValAbs_opt = par_min+0.5*diag(DELTA_opt+ones(num_uncertain_par,1))*D;
    disp('-----------------------------------------------------------------------------------------');
    disp('|  Name            |  Initial Value  |  Optimal Value  |  Lower Bound  |  Higher Bound  |');
    disp('-----------------------------------------------------------------------------------------');
    for rr = 1:num_uncertain_par
        row = ['|  ',ext_lftfun.DeltaSym{1+rr,1}];
        for cc=1:16-length(ext_lftfun.DeltaSym{1+rr,1})
            row = [row,32, ''];
        end;
        row = [row, '|  ', num2str(par_ValAbs_0(rr))];
        for cc=1:15-length(num2str(par_ValAbs_0(rr)))
            row = [row,32, ''];
        end;
        row = [row, '|  ', num2str(par_ValAbs_opt(rr))];
        for cc=1:15-length(num2str(par_ValAbs_opt(rr)))
            row = [row,32, ''];
        end;
        row = [row, '|  ', num2str(delta_min(rr))];
        for cc=1:13-length(num2str(delta_min(rr)))
            row = [row,32, ''];
        end;
        row = [row, '|  ', num2str(delta_max(rr))];
        for cc=1:14-length(num2str(delta_max(rr)))
            row = [row,32, ''];
        end;        
        row = [row, '|'];
        disp(row);
    end;
    disp('-----------------------------------------------------------------------------------------');
