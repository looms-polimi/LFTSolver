function LFTCreator(tag,...
                    A_ode_sim,...
                    B_csi_ode_sim,...
                    B_u_ode_sim,...
                    C_ode_sim,...
                    D_ode_sim,...
                    Theta_ode_sim,...
                    p_sim,...
                    u_sim,...
                    x_sim,...
                    p_min,...
                    p_max,...
                    u_names,...
                    x_names,...
                    y_names,...
                    nominal,...
                    A_temp_sim,...
                    B_u_temp_sim,...
                    pnl_sim,...
                    z_sim,...
                    pnl_min,...
                    pnl_max)
for i=1:length(p_sim)
    eval(['syms ',char(p_sim(i))]);
end;
for i=1:length(pnl_sim)
    eval(['syms ',char(pnl_sim(i))]);
end;
for i=1:length(u_sim)
    eval(['syms ',char(u_sim(i))]);
end;
for i=1:length(x_sim)
    eval(['syms ',char(x_sim(i))]);
end;
for i=1:length(z_sim)
    eval(['syms ',char(z_sim(i))]);
end;
nx = size(A_ode_sim,1);
nc = size(B_csi_ode_sim,2);
nu = size(B_u_ode_sim,2);
ny = size(C_ode_sim,1);
n_parameters = length(p_sim);
nzwnl = length(pnl_sim);
DeltaSym = cat(1,{'parName', 'indStartDiag', 'indStopDiag', 'LowerBound', 'HigherBound', 'toIdentify', 'lb', 'ub'},cell(n_parameters,8));
tag = [tag,'_LFT'];


%% costruzione elenchi parametri per matrici ode
A_ode_par = cell(nx,nx);
for i=1:nx
    for j=1:nx
        A_ode_par{i,j} = symvar(A_ode_sim(i,j));   
    end;
end;
B_csi_ode_par = cell(nx,nc);
for i=1:nx
    for j=1:nc
        B_csi_ode_par{i,j} = symvar(B_csi_ode_sim(i,j));   
    end;
end;
B_u_ode_par = cell(nx,nu);
for i=1:nx
    for j=1:nu
        B_u_ode_par{i,j} = symvar(B_u_ode_sim(i,j));   
    end;
end;
C_ode_par = cell(ny,nx);
for i=1:ny
    for j=1:nx
        C_ode_par{i,j} = symvar(C_ode_sim(i,j));   
    end;
end;
D_ode_par = cell(ny,nu);
for i=1:ny
    for j=1:nu
        D_ode_par{i,j} = symvar(D_ode_sim(i,j));   
    end;
end;


%% Costruzione elenchi variabili in equazioni non lineari
Theta_ode_par = cell(nc,1);
for i=1:nc
    Theta_ode_par{i} = symvar(Theta_ode_sim(i));   
end;


%% Costruzione matrice A B2 B3 C3
A = A_ode_sim;
B2 = B_csi_ode_sim;
B3 = B_u_ode_sim;
C3 = C_ode_sim;
D33 = D_ode_sim;


%% Costruzione matrici C1 D12 D13
C1 = sym([]);
D12 = sym([]);
D13 = sym([]);
nzw = 0;  % indice parametro in matrice DELTA
state_eq_index = [];
output_eq_index = [];
macroTerm_coord = [];
for p=1:n_parameters
    par_name = ['p',num2str(p)];
    par_sym = sym(par_name);
    delta_par_sym = sym(['p',num2str(p),'_delta']);
    for i=1:nx
        for j=1:nx
            if ~isempty(A_ode_par{i,j})
                term_position = find(A_ode_par{i,j}==par_sym);
                if ~isempty(term_position)
                    C1 = cat(1,C1,sym(zeros(1,nx)));
                    D12 = cat(1,D12,sym(zeros(1,nc)));
                    D13 = cat(1,D13,sym(zeros(1,nu)));
                    if term_position==1
                        expr = char(A_ode_sim(i,j));
                        expr = subs(expr, par_sym, delta_par_sym);
                        nzw = nzw+1;
                        state_eq_index = cat(1,state_eq_index,i);
                        output_eq_index = cat(1,output_eq_index,NaN);
                        C1(nzw,j) = expr;
                    else % da estendere a casi più complessi che per ora non servono (controllo segni e parentesi)
                        nzw = nzw+1;
                        expr_new = char(A_ode_sim(i,j));
                        ind_new = strfind(expr_new,char(par_sym));
                        if strcmp(expr_new(ind_new(1)-1),'*')
                            ff = 1;
                            factor_str = char();
                            while ~isnan(str2double(expr_new(ind_new(1)-1-ff)))
                                factor_str = [expr_new(ind_new(1)-1-ff),factor_str];
                                ff = ff+1;
                            end;
                            if strcmp(expr_new(ind_new(1)-1-ff),'(')||strcmp(expr_new(ind_new(1)-1-ff),' ')
                                factor = str2double(factor_str);
                            else
                                factor = 1;
                            end;
                        else
                            factor = 1;
                        end;
                        state_eq_index = cat(1,state_eq_index,i);
                        output_eq_index = cat(1,output_eq_index,NaN);
                        C1(nzw,j) = sym(factor)*delta_par_sym;
                    end;
                    macroTerm_coord = cat(1,macroTerm_coord,[i,j,1]);
                    if isempty(DeltaSym{1+p,1})
                        DeltaSym{1+p,1} = par_name;
                        if p==1
                            DeltaSym{1+p,2}=1;
                            DeltaSym{1+p,3}=DeltaSym{1+p,2};
                        else
                            DeltaSym{1+p,2}=DeltaSym{1+p-1,3}+1;
                            DeltaSym{1+p,3}=DeltaSym{1+p,2};
                        end;
                    else
                        DeltaSym{1+p,3} = DeltaSym{1+p,3}+1;
                    end;
                end;
            end;
        end
    end
    for i=1:nx
        for j=1:nc
            if ~isempty(B_csi_ode_par{i,j})
                term_position = find(B_csi_ode_par{i,j}==par_sym);
                if ~isempty(term_position)
                    C1 = cat(1,C1,sym(zeros(1,nx)));
                    D12 = cat(1,D12,sym(zeros(1,nc)));
                    D13 = cat(1,D13,sym(zeros(1,nu)));
                    if term_position==1
                        expr = char(B_csi_ode_sim(i,j));       
                        expr = subs(expr, par_sym, delta_par_sym);
                        nzw = nzw+1;
                        state_eq_index = cat(1,state_eq_index,i);
                        output_eq_index = cat(1,output_eq_index,NaN);
                        D12(nzw,j) = expr;
                    else % da estendere a casi più complessi che per ora non servono (controllo segni e parentesi)
                        nzw = nzw+1;
                        expr_new = char(B_csi_ode_sim(i,j));
                        ind_new = strfind(expr_new,char(par_sym));
                        if strcmp(expr_new(ind_new(1)-1),'*')
                            ff = 1;
                            factor_str = char();
                            while ~isnan(str2double(expr_new(ind_new(1)-1-ff)))
                                factor_str = [expr_new(ind_new(1)-1-ff),factor_str];
                                ff = ff+1;
                            end;
                            if strcmp(expr_new(ind_new(1)-1-ff),'(')||strcmp(expr_new(ind_new(1)-1-ff),' ')
                                factor = str2double(factor_str);
                            else
                                factor = 1;
                            end;
                        else
                            factor = 1;
                        end;
                        state_eq_index = cat(1,state_eq_index,i);
                        output_eq_index = cat(1,output_eq_index,NaN);
                        D12(nzw,j) = sym(factor)*delta_par_sym;
                    end;
                    macroTerm_coord = cat(1,macroTerm_coord,[i,j,2]);
                    if isempty(DeltaSym{1+p,1})
                        DeltaSym{1+p,1} = par_name;
                        if p==1
                            DeltaSym{1+p,2}=1;
                            DeltaSym{1+p,3}=DeltaSym{1+p,2};
                        else
                            DeltaSym{1+p,2}=DeltaSym{1+p-1,3}+1;
                            DeltaSym{1+p,3}=DeltaSym{1+p,2};
                        end;
                    else
                        DeltaSym{1+p,3} = DeltaSym{1+p,3}+1;
                    end;
                end;
            end;
        end
    end
    for i=1:nx
        for j=1:nu
            if ~isempty(B_u_ode_par{i,j})
                term_position = find(B_u_ode_par{i,j}==par_sym);
                if ~isempty(term_position)
                    C1 = cat(1,C1,sym(zeros(1,nx)));
                    D12 = cat(1,D12,sym(zeros(1,nc)));
                    D13 = cat(1,D13,sym(zeros(1,nu)));
                    if term_position==1
                        expr = char(B_u_ode_sim(i,j));
                        expr = subs(expr, par_sym, delta_par_sym);
                        nzw = nzw+1;
                        state_eq_index = cat(1,state_eq_index,i);
                        output_eq_index = cat(1,output_eq_index,NaN);
                        D13(nzw,j) = expr;
                    else % da estendere a casi più complessi che per ora non servono (controllo segni e parentesi)
                        nzw = nzw+1;
                        expr_new = char(B_u_ode_sim(i,j));
                        ind_new = strfind(expr_new,char(par_sym));
                        if strcmp(expr_new(ind_new(1)-1),'*')
                            ff = 1;
                            factor_str = char();
                            while ~isnan(str2double(expr_new(ind_new(1)-1-ff)))
                                factor_str = [expr_new(ind_new(1)-1-ff),factor_str];
                                ff = ff+1;
                            end;
                            if strcmp(expr_new(ind_new(1)-1-ff),'(')||strcmp(expr_new(ind_new(1)-1-ff),' ')
                                factor = str2double(factor_str);
                            else
                                factor = 1;
                            end;
                        else
                            factor = 1;
                        end;
                        state_eq_index = cat(1,state_eq_index,i);
                        output_eq_index = cat(1,output_eq_index,NaN);
                        D13(nzw,j) = sym(factor)*delta_par_sym;
                    end;
                    macroTerm_coord = cat(1,macroTerm_coord,[i,j,3]);
                    if isempty(DeltaSym{1+p,1})
                        DeltaSym{1+p,1} = par_name;
                        if p==1
                            DeltaSym{1+p,2}=1;
                            DeltaSym{1+p,3}=DeltaSym{1+p,2};
                        else
                            DeltaSym{1+p,2}=DeltaSym{1+p-1,3}+1;
                            DeltaSym{1+p,3}=DeltaSym{1+p,2};
                        end;
                    else
                        DeltaSym{1+p,3} = DeltaSym{1+p,3}+1;
                    end;
                end;
            end;
        end
    end
    for i=1:ny
        for j=1:nx
            if ~isempty(C_ode_par{i,j})
                term_position = find(C_ode_par{i,j}==par_sym);
                if ~isempty(term_position)
                    C1 = cat(1,C1,sym(zeros(1,nx)));
                    D12 = cat(1,D12,sym(zeros(1,nc)));
                    D13 = cat(1,D13,sym(zeros(1,nu)));
                    if term_position==1
                        expr = char(C_ode_sim(i,j));
                        expr = subs(expr, par_sym, delta_par_sym);
                        nzw = nzw+1;
                        output_eq_index = cat(1,output_eq_index,i);
                        state_eq_index = cat(1,state_eq_index,NaN);
                        C1(nzw,j) = expr;
                    else % da estendere a casi più complessi che per ora non servono (controllo segni e parentesi)
                        nzw = nzw+1;
                        expr_new = char(C_ode_sim(i,j));
                        ind_new = strfind(expr_new,char(par_sym));
                        if strcmp(expr_new(ind_new(1)-1),'*')
                            ff = 1;
                            factor_str = char();
                            while ~isnan(str2double(expr_new(ind_new(1)-1-ff)))
                                factor_str = [expr_new(ind_new(1)-1-ff),factor_str];
                                ff = ff+1;
                            end;
                            if strcmp(expr_new(ind_new(1)-1-ff),'(')||strcmp(expr_new(ind_new(1)-1-ff),' ')
                                factor = str2double(factor_str);
                            else
                                factor = 1;
                            end;
                        else
                            factor = 1;
                        end;
                        output_eq_index = cat(1,output_eq_index,i);
                        state_eq_index = cat(1,state_eq_index,NaN);
                        C1(nzw,j) = sym(factor)*delta_par_sym;
                    end;
                    macroTerm_coord = cat(1,macroTerm_coord,[i,j,4]);
                    if isempty(DeltaSym{1+p,1})
                        DeltaSym{1+p,1} = par_name;
                        if p==1
                            DeltaSym{1+p,2}=1;
                            DeltaSym{1+p,3}=DeltaSym{1+p,2};
                        else
                            DeltaSym{1+p,2}=DeltaSym{1+p-1,3}+1;
                            DeltaSym{1+p,3}=DeltaSym{1+p,2};
                        end;
                    else
                        DeltaSym{1+p,3} = DeltaSym{1+p,3}+1;
                    end;
                end;
            end;
        end
    end
    for i=1:ny
        for j=1:nu
            if ~isempty(D_ode_par{i,j})
                term_position = find(D_ode_par{i,j}==par_sym);
                if ~isempty(term_position)
                    C1 = cat(1,C1,sym(zeros(1,nx)));
                    D12 = cat(1,D12,sym(zeros(1,nc)));
                    D13 = cat(1,D13,sym(zeros(1,nu)));
                    if term_position==1
                        expr = char(D_ode_sim(i,j));
                        expr = subs(expr, par_sym, delta_par_sym);
                        nzw = nzw+1;
                        output_eq_index = cat(1,output_eq_index,i);
                        state_eq_index = cat(1,state_eq_index,NaN);
                        C1(nzw,j) = expr;
                    else % da estendere a casi più complessi che per ora non servono (controllo segni e parentesi)
                        nzw = nzw+1;
                        expr_new = char(D_ode_sim(i,j));
                        ind_new = strfind(expr_new,char(par_sym));
                        if strcmp(expr_new(ind_new(1)-1),'*')
                            ff = 1;
                            factor_str = char();
                            while ~isnan(str2double(expr_new(ind_new(1)-1-ff)))
                                factor_str = [expr_new(ind_new(1)-1-ff),factor_str];
                                ff = ff+1;
                            end;
                            if strcmp(expr_new(ind_new(1)-1-ff),'(')||strcmp(expr_new(ind_new(1)-1-ff),' ')
                                factor = str2double(factor_str);
                            else
                                factor = 1;
                            end;
                        else
                            factor = 1;
                        end;
                        output_eq_index = cat(1,output_eq_index,i);
                        state_eq_index = cat(1,state_eq_index,NaN);
                        C1(nzw,j) = sym(factor)*delta_par_sym;
                    end;
                    macroTerm_coord = cat(1,macroTerm_coord,[i,j,4]);
                    if isempty(DeltaSym{1+p,1})
                        DeltaSym{1+p,1} = par_name;
                        if p==1
                            DeltaSym{1+p,2}=1;
                            DeltaSym{1+p,3}=DeltaSym{1+p,2};
                        else
                            DeltaSym{1+p,2}=DeltaSym{1+p-1,3}+1;
                            DeltaSym{1+p,3}=DeltaSym{1+p,2};
                        end;
                    else
                        DeltaSym{1+p,3} = DeltaSym{1+p,3}+1;
                    end;
                end;
            end;
        end
    end

end;


%% Costruzione matrici C2 D21 D23
C2 = sym([]);
D21 = sym([]);
D23 = sym([]);
equation = sym([]);
no= 0;   % indice riga per ingresso non lineare omega
i_var_temp = 0;
for i=1:nc
    terms = Theta_ode_par{i};
    equation = cat(1,equation,Theta_ode_sim(i));
    for j=1:length(terms)
        term_sym = terms(j);
        term_name = char(term_sym);
        term_type = term_name(1);                    % estendere anche all'ingresso di eventuali w
        term_ind = str2double(term_name(2:end));
        C2 = cat(1,C2,sym(zeros(1,nx)));
        D21 = cat(1,D21,sym(zeros(1,nzw+nzwnl)));
        D23 = cat(1,D23,sym(zeros(1,nu)));
        no = no+1;
        switch term_type                             % estendere anche all'ingresso di eventuali w
            case 'x'
                C2(no,term_ind) = 1;
            case 'u'
                D23(no,term_ind) = 1;
            case 'z'
                i_var_temp = i_var_temp + 1;
                if ~all(A_temp_sim(i_var_temp,:)==sym(zeros(1,nx)))
                    ind_col_pnl = find(~(A_temp_sim(i_var_temp,:)==0));
                    name_pnl = char(A_temp_sim(i_var_temp,ind_col_pnl));
                    vector_pnl = A_temp_sim(i_var_temp,:);
                    vector_pnl = subs(vector_pnl, A_temp_sim(i_var_temp,ind_col_pnl), sym([name_pnl,'_delta']));
                    C1 = cat(1,C1,vector_pnl);
                    D12 = cat(1,D12,sym(zeros(1,nc)));
                    D13 = cat(1,D13,sym(zeros(1,nu)));
                    C2(no,:) = A_temp_sim(i_var_temp,:);
                    D21(no,nzw+i_var_temp) = 1;   
                elseif ~all(B_u_temp_sim(i_var_temp,:)==sym(zeros(1,nu)))
                    ind_col_pnl = find(~(B_u_temp_sim(i_var_temp,:)==0));
                    name_pnl = char(B_u_temp_sim(i_var_temp,ind_col_pnl));
                    vector_pnl = B_u_temp_sim(i_var_temp,:);
                    vector_pnl = subs(vector_pnl, B_u_temp_sim(i_var_temp,ind_col_pnl), sym([name_pnl,'_delta']));
                    C1 = cat(1,C1,sym(zeros(1,nx)));
                    D12 = cat(1,D12,sym(zeros(1,nc)));
                    D13 = cat(1,D13,vector_pnl);
                    D21(no,nzw+i_var_temp) = 1;
                    D23(no,:) = B_u_temp_sim(i_var_temp,:);
                else
                    error('wrong definition for temporary matrices in specifications');
                end; 
            otherwise
                error('bad value for non-linear input arguments');
        end;
        omega_sym = sym(['omega_',int2str(no)]);
        equation(i) = subs(equation(i), term_sym, omega_sym);
    end;
end;


%%  costruzione Theta
nomOut1 = [tag, '_Theta.m'];
g1 = fopen(nomOut1, 'wt');
fprintf(g1, '%s\n', ['function  y = ', tag,'_Theta(omega)']);
for i=1:no
     fprintf(g1, '%s\n', ['omega_',int2str(i),'=omega(',int2str(i),');']);
end
fprintf(g1, '%s\n', 'y=[');
for i=1:nc
     fprintf(g1, '%s\n', char(equation(i)));
end
fprintf(g1, '%s\n', '];');
fclose(g1);


%%  costruzione dThetadOmega
omega = [];
for i=1:no
     eval(['syms omega_',int2str(i)]);
     omega = cat(1,omega,eval(['omega_',int2str(i)]));
end
y = equation;
dydo = jacobian(y, omega);
nomOut2 = [tag, '_dThetadOmega.m'];
g2 = fopen(nomOut2, 'wt');
fprintf(g2, '%s\n', ['function  y = ', tag,'_dThetadOmega(omega)']);
for i=1:no
     fprintf(g2, '%s\n', ['omega_',int2str(i),'=omega(',int2str(i),');']);
end
fprintf(g2, '%s\n', 'y=[');
for i=1:nc
     dydo_row_i_test = char(dydo(i,:));
     if length(dydo_row_i_test)>1
         fprintf(g2, '%s\n', dydo_row_i_test(10:end-3));
     else
         fprintf(g2, '%s\n', dydo_row_i_test(1:end));
     end;
end
fprintf(g2, '%s\n', '];');
fclose(g2);


%% Costruzione matrici D22 D32 C3
D22 = sym(zeros(no,nc));
D32 = sym(zeros(ny,nc));


%% Costruzione matrici B1 D11 D31
B1 = sym(zeros(nx,nzw));
D11 = sym(zeros(nzw,nzw));
D31 = sym(zeros(ny,nzw+nzwnl));
for i=1:nzw
    for j=1:nx
        par_expr = C1(i,j);
        terms = symvar(par_expr);
        if length(terms)>1                  % non è un singolo parametro allora entra in shell
            ind_par_in_shell = [];
            for k=1:length(terms)           % per ora è così, poi questo andrà rimpiazzato con un ciclo che studia le parentesi e i segni
                par_sym = terms(k);
                par_name = char(par_sym);
                ind_ = strfind(par_name,'_delta');
                if isempty(ind_)
                    ind_par_in_shell = cat(1,ind_par_in_shell,str2double(par_name(2:end)));
                else
                    par_delta_out_shell = par_sym;
                    sign_par_out_shell = +1;
                    par_expr_string = char(par_expr);
                    if strcmp(par_expr_string(1),'-')
                        sign_par_out_shell = -1;
                    end;
                end;
                par_expr_string = char(par_expr);
                if strcmp(par_expr_string(1),'-')
                    ff = 2;
                else
                    ff = 1;
                end;
                factor = 1;
                factor_str = char();
                while ~isnan(str2double(par_expr_string(ff)))
                    factor_str = [factor_str,par_expr_string(ff)];
                    factor = str2double(factor_str);
                    ff =ff+1;
                end;
            end;
            ind_col = [];
            for cc=1:length(ind_par_in_shell)
                for ii=1:nzw
                    par_expr = C1(ii,j);
                    terms = symvar(par_expr);
                    if length(terms)==1
                        par_name_delta_shell = ['p',num2str(ind_par_in_shell(cc)),'_delta'];
                        if ~isempty(find(terms==sym(par_name_delta_shell)))&&all(macroTerm_coord(i,:)==macroTerm_coord(ii,:))
                            ind_col = ii;
                            break;
                        end;
                    end;
                end;
                par_delta_out_shell_string = char(par_delta_out_shell);
                par_out_shell = sym(par_delta_out_shell_string(1:ind_-1));
                if (macroTerm_coord(i,3)==1) %||(macroTerm_coord(i,3)==2)||(macroTerm_coord(i,3)==3) per gli altri casi (Y=f(csi,u))
                    D11(i,ind_col) = D11(i,ind_col) + sym(factor)*sign_par_out_shell*par_delta_out_shell;
                    B1(state_eq_index(i),ind_col) = B1(state_eq_index(i),ind_col) + sym(factor)*sign_par_out_shell*par_out_shell;
                    B1(state_eq_index(i),i) = 1;
                elseif (macroTerm_coord(i,3)==4) %||(macroTerm_coord(i,3)==5)||(macroTerm_coord(i,3)==6) per gli altri casi (Y=f(csi,u))
                    D11(i,ind_col) = D11(i,ind_col) + sym(factor)*sign_par_out_shell*par_delta_out_shell;
                    D31(output_eq_index(i),ind_col) = D31(output_eq_index(i),ind_col) + sym(factor)*sign_par_out_shell*par_out_shell;
                    D31(output_eq_index(i),i) = 1;
                else
                    error('wrong value for macroTerm_coord(#,#,v)');
                end;
            end;
        else
            if length(terms)==1
                if (macroTerm_coord(i,3)==1) %||(macroTerm_coord(i,3)==2)||(macroTerm_coord(i,3)==3) per gli altri casi (Y=f(csi,u))
                    term_list = A_ode_par{macroTerm_coord(i,1),macroTerm_coord(i,2)};
                    if length(term_list)==1
                        B1(state_eq_index(i),i) = 1;
                    end;
                elseif (macroTerm_coord(i,3)==4) %||(macroTerm_coord(i,3)==5)||(macroTerm_coord(i,3)==6) per gli altri casi (Y=f(csi,u))
                    term_list = C_ode_par{macroTerm_coord(i,1),macroTerm_coord(i,2)};
                    if length(term_list)==1
                        D31(output_eq_index(i),i) = 1;
                    end;
                else
                    error('wrong value for macroTerm_coord(#,#,v)');
                end;    
            end;
        end;
    end;
end;
for i=1:nzw
    for j=1:nc
        par_expr = D12(i,j);
        terms = symvar(par_expr);
        if length(terms)>1                  % non è un singolo parametro allora entra in shell
            ind_par_in_shell = [];
            for k=1:length(terms)           % per ora è così, poi questo andrà rimpiazzato con un ciclo che studia le parentesi e i segni
                par_sym = terms(k);
                par_name = char(par_sym);
                ind_ = strfind(par_name,'_delta');
                if isempty(ind_)
                    ind_par_in_shell = cat(1,ind_par_in_shell,str2double(par_name(2:end)));
                else
                    par_delta_out_shell = par_sym;
                    sign_par_out_shell = +1;
                    par_expr_string = char(par_expr);
                    if strcmp(par_expr_string(1),'-')
                        sign_par_out_shell = -1;
                    end;
                end;
                par_expr_string = char(par_expr);
                if strcmp(par_expr_string(1),'-')
                    ff = 2;
                else
                    ff = 1;
                end;
                factor = 1;
                factor_str = char();
                while ~isnan(str2double(par_expr_string(ff)))
                    factor_str = [factor_str,par_expr_string(ff)];
                    factor = str2double(factor_str);
                    ff =ff+1;
                end; 
            end;
            ind_col = [];
            for cc=1:length(ind_par_in_shell)
                for ii=1:nzw
                    par_expr = D12(ii,j);
                    terms = symvar(par_expr);
                    if length(terms)==1
                        par_name_delta_shell = ['p',num2str(ind_par_in_shell(cc)),'_delta'];
                        if ~isempty(find(terms==sym(par_name_delta_shell)))&&all(macroTerm_coord(i,:)==macroTerm_coord(ii,:))
                            ind_col = ii;
                            break;
                        end;
                    end;
                end;
                par_delta_out_shell_string = char(par_delta_out_shell);
                par_out_shell = sym(par_delta_out_shell_string(1:ind_-1));
                D11(i,ind_col) = D11(i,ind_col) + sym(factor)*sign_par_out_shell*par_delta_out_shell;
                B1(state_eq_index(i),ind_col) = B1(state_eq_index(i),ind_col) + sym(factor)*sign_par_out_shell*par_out_shell;
                B1(state_eq_index(i),i) = 1;
            end;
        else
            if length(terms)==1
                term_list = B_csi_ode_par{macroTerm_coord(i,1),macroTerm_coord(i,2)};
                if length(term_list)==1
                    B1(state_eq_index(i),i) = 1;
                end;
            end;
        end;
    end;
end;
for i=1:nzw
    for j=1:nu
        par_expr = D13(i,j);
        terms = symvar(par_expr);
        if length(terms)>1                  % non è un singolo parametro allora entra in shell
            ind_par_in_shell = [];
            for k=1:length(terms)           % per ora è così, poi questo andrà rimpiazzato con un ciclo che studia le parentesi e i segni
                par_sym = terms(k);
                par_name = char(par_sym);
                ind_ = strfind(par_name,'_delta');
                if isempty(ind_)
                    ind_par_in_shell = cat(1,ind_par_in_shell,str2double(par_name(2:end)));
                else
                    par_delta_out_shell = par_sym;
                    sign_par_out_shell = +1;
                    par_expr_string = char(par_expr);
                    if strcmp(par_expr_string(1),'-')
                        sign_par_out_shell = -1;
                    end;
                end;
                par_expr_string = char(par_expr);
                if strcmp(par_expr_string(1),'-')
                    ff = 2;
                else
                    ff = 1;
                end;
                factor = 1;
                factor_str = char();
                while ~isnan(str2double(par_expr_string(ff)))
                    factor_str = [factor_str,par_expr_string(ff)];
                    factor = str2double(factor_str);
                    ff =ff+1;
                end; 
            end;
            ind_col = [];
            for cc=1:length(ind_par_in_shell)
                for ii=1:nzw
                    par_expr = D13(ii,j);
                    terms = symvar(par_expr);
                    if length(terms)==1
                        par_name_delta_shell = ['p',num2str(ind_par_in_shell(cc)),'_delta'];
                        if ~isempty(find(terms==sym(par_name_delta_shell)))&&all(macroTerm_coord(i,:)==macroTerm_coord(ii,:))
                            ind_col = ii;
                            break;
                        end;
                    end;
                end;
                par_delta_out_shell_string = char(par_delta_out_shell);
                par_out_shell = sym(par_delta_out_shell_string(1:ind_-1));
                D11(i,ind_col) = D11(i,ind_col) + sym(factor)*sign_par_out_shell*par_delta_out_shell;
                B1(state_eq_index(i),ind_col) = B1(state_eq_index(i),ind_col) + sym(factor)*sign_par_out_shell*par_out_shell;
                B1(state_eq_index(i),i) = 1;
            end;
        else
            if length(terms)==1
                term_list = B_u_ode_par{macroTerm_coord(i,1),macroTerm_coord(i,2)};
                if length(term_list)==1
                    B1(state_eq_index(i),i) = 1;
                end;
            end;
        end;
    end;
end;


%% contornamento matrici con zeri
B1 = cat(2, B1, sym(zeros(nx,nzwnl)));
D11 = cat(2, D11, sym(zeros(nzw,nzwnl)));
D11 = cat(1, D11, sym(zeros(nzwnl,nzw+nzwnl)));


%% valutazione e assemblaggio struttura forma LTI, completamento tabelle e informazioni di sparsità
p_0 = zeros(n_parameters,1);
p_delta = zeros(n_parameters,1);
for pp=1:n_parameters
    DeltaSym{1+pp,4} = p_min(pp);
    DeltaSym{1+pp,5} = p_max(pp);
    DeltaSym{1+pp,6} = 1;
    p_0(pp) = 0.5*(p_max(pp)+p_min(pp));
    p_delta(pp) = 0.5*(p_max(pp)-p_min(pp));
    eval(['p',num2str(pp),' = p_0(',num2str(pp),');']);
    eval(['p',num2str(pp),'_delta = p_delta(',num2str(pp),');']);
end;
pnl_0 = zeros(nzwnl,1);
pnl_delta = zeros(nzwnl,1);
for pp=1:nzwnl
    DeltaSym = cat(1,DeltaSym,cell(1,8));
    DeltaSym{1+n_parameters+pp,1} = ['pnl',num2str(pp)];
    DeltaSym{1+n_parameters+pp,2} = DeltaSym{1+n_parameters,3}+pp;
    DeltaSym{1+n_parameters+pp,3} = DeltaSym{1+n_parameters,3}+pp;
    DeltaSym{1+n_parameters+pp,4} = pnl_min(pp);
    DeltaSym{1+n_parameters+pp,5} = pnl_max(pp);
    DeltaSym{1+n_parameters+pp,6} = 1;
    pnl_0(pp) = 0.5*(pnl_max(pp)+pnl_min(pp));
    pnl_delta(pp) = 0.5*(pnl_max(pp)-pnl_min(pp));
    eval(['pnl',num2str(pp),' = pnl_0(',num2str(pp),');']);
    eval(['pnl',num2str(pp),'_delta = pnl_delta(',num2str(pp),');']);
end;
DeltaVal = zeros(nzw+nzwnl,nzw+nzwnl);
A_num= eval(A);
B1_num= eval(B1);
B2_num= eval(B2);
B3_num= eval(B3);
C1_num= eval(C1);
C2_num= eval(C2);
C3_num= eval(C3);
D11_num= eval(D11);
D12_num= eval(D12);
D13_num= eval(D13);
D21_num= eval(D21);
D22_num= eval(D22);
D23_num= eval(D23);
D31_num= eval(D31);
D32_num= eval(D32);
D33_num= eval(D33);
LTI = struct('A',A_num,     'B1',B1_num,    'B2',B2_num,    'B3',B3_num,...
             'C1',C1_num,   'D11',D11_num,  'D12',D12_num,  'D13',D13_num,...
             'C2',C2_num,   'D21',D21_num,  'D22',D22_num,  'D23',D23_num,...
             'C3',C3_num,   'D31',D31_num,  'D32',D32_num,  'D33',D33_num);
A_Sparsity = sparsityMatrix(A_num);
B1_Sparsity = sparsityMatrix(B1_num);
B2_Sparsity = sparsityMatrix(B2_num);
B3_Sparsity = sparsityMatrix(B3_num);
C1_Sparsity = sparsityMatrix(C1_num);
C2_Sparsity = sparsityMatrix(C2_num);
C3_Sparsity = sparsityMatrix(C3_num);
D11_Sparsity = sparsityMatrix(D11_num);
D12_Sparsity = sparsityMatrix(D12_num);
D13_Sparsity = sparsityMatrix(D13_num);
D21_Sparsity = sparsityMatrix(D21_num);
D22_Sparsity = sparsityMatrix(D22_num);
D23_Sparsity = sparsityMatrix(D23_num);
D31_Sparsity = sparsityMatrix(D31_num);
D32_Sparsity = sparsityMatrix(D32_num);
D33_Sparsity = sparsityMatrix(D33_num);
dThetadOmega_Sparsity = true(size(dydo));

[row_zero_elements, col_zero_elements] = find(dydo==0);
if ~isempty(find(dydo==0))
    for ind_coord = 1:size(row_zero_elements,1)
        dThetadOmega_Sparsity(row_zero_elements(ind_coord),col_zero_elements(ind_coord))=false;
    end;
end;

VCT8_Sparsity = rowColProductBoolean(invBoolean([((diag(true(size(C1_num,1),1)))|D11_Sparsity), rowColProductBoolean(D12_Sparsity,dThetadOmega_Sparsity)
                                                 D21_Sparsity,                              ((diag(true(size(C2_num,1),1)))|rowColProductBoolean(D22_Sparsity,dThetadOmega_Sparsity))]),[C1_Sparsity, D13_Sparsity, D11_Sparsity
                                                                                                                                                                                     C2_Sparsity, D23_Sparsity, D21_Sparsity]);                                                 
VCT9_Sparsity = ([A_Sparsity,  B1_Sparsity
                  C3_Sparsity, D31_Sparsity]|rowColProductBoolean([B1_Sparsity,  rowColProductBoolean(B2_Sparsity,dThetadOmega_Sparsity)
                                                                   D31_Sparsity, rowColProductBoolean(D32_Sparsity,dThetadOmega_Sparsity)],[VCT8_Sparsity(:,1:size(A_num,2)),VCT8_Sparsity(:,size(A_num,2)+size(D13_num,2)+1:end)]));

VCT10_Sparsity = ([A_Sparsity, B3_Sparsity
                  C3_Sparsity, D33_Sparsity]|rowColProductBoolean([B1_Sparsity,  rowColProductBoolean(B2_Sparsity,dThetadOmega_Sparsity)
                                                                   D31_Sparsity, rowColProductBoolean(D32_Sparsity,dThetadOmega_Sparsity)],VCT8_Sparsity(:,1:size(A_num,2)+size(D13_num,2))));
                                      
SparsityInformations = struct();

[row,col] = find(VCT9_Sparsity(1:size(A_num,1),1:size(A_num,2)));
SparsityInformations.VCT9{1,1}.NumNonZeroElements = length(row);
SparsityInformations.VCT9{1,1}.row = row;
SparsityInformations.VCT9{1,1}.col = col;

[row,col] = find(VCT9_Sparsity(size(A_num,1)+1:end,1:size(A_num,2)));
SparsityInformations.VCT9{2,1}.NumNonZeroElements = length(row);
SparsityInformations.VCT9{2,1}.row = row;
SparsityInformations.VCT9{2,1}.col = col;

for i=1:n_parameters+nzwnl
    ind_start_z = DeltaSym{1+i,2};
    ind_stop_z = DeltaSym{1+i,3};

    [row,col] = find(VCT9_Sparsity(1:size(A_num,1),size(A_num,2)+ind_start_z:size(A_num,2)+ind_stop_z));
    SparsityInformations.VCT9{1,1+i}.NumNonZeroElements = length(row);
    SparsityInformations.VCT9{1,1+i}.row = row;
    SparsityInformations.VCT9{1,1+i}.col = col;
    
    [row,col] = find(VCT9_Sparsity(size(A_num,1)+1:end,size(A_num,2)+ind_start_z:size(A_num,2)+ind_stop_z));
    SparsityInformations.VCT9{2,1+i}.NumNonZeroElements = length(row);
    SparsityInformations.VCT9{2,1+i}.row = row;
    SparsityInformations.VCT9{2,1+i}.col = col;
end;

[row,col] = find(VCT10_Sparsity(1:size(A_num,1),1:size(A_num,2)));
SparsityInformations.VCT10{1,1}.NumNonZeroElements = length(row);
SparsityInformations.VCT10{1,1}.row = row;
SparsityInformations.VCT10{1,1}.col = col;

[row,col] = find(VCT10_Sparsity(size(A_num,1)+1:end,1:size(A_num,2)));
SparsityInformations.VCT10{2,1}.NumNonZeroElements = length(row);
SparsityInformations.VCT10{2,1}.row = row;
SparsityInformations.VCT10{2,1}.col = col;

for i=1:size(B3_num,2)
    [row,col] = find(VCT10_Sparsity(1:size(A_num,1),size(A_num,2)+i));
    SparsityInformations.VCT10{1,1+i}.NumNonZeroElements = length(row);
    SparsityInformations.VCT10{1,1+i}.row = row;
    SparsityInformations.VCT10{1,1+i}.col = col;

    [row,col] = find(VCT10_Sparsity(size(A_num,1)+1:end,size(A_num,2)+i));
    SparsityInformations.VCT10{2,1+i}.NumNonZeroElements = length(row);
    SparsityInformations.VCT10{2,1+i}.row = row;
    SparsityInformations.VCT10{2,1+i}.col = col;
end;

for i=1:n_parameters+nzwnl
    DeltaSym{1+i,7} = -1.9;
    DeltaSym{1+i,8} = Inf;
end;


%% struttura dati finale
lftfun = struct('LTI', LTI,...
                'DeltaSym', {DeltaSym},...
                'DeltaVal', DeltaVal,...
                'Theta', eval(['@', tag, '_Theta']),...
                'dThetadOmega', eval(['@', tag, '_dThetadOmega']),...
                'Sparsity', SparsityInformations,...
                'ISO', struct('InputVariablesNames', {u_names},...
                              'StateVariablesNames', {x_names},...
                              'OutputVariablesNames', {y_names},...
                              'Nominal', nominal));                    
eval([tag, ' = lftfun']);
save([tag, '.mat'], tag); 



end



