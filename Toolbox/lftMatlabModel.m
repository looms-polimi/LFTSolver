function lftMatlabModel(nome_m,nome_xml)
%--------------------------------------------------------------------------
% Generazione automatica della struttura LFT per LFTsolver a partire dal
% file .m generato da ECLIPSE e dal file .xml generato da OMC
% Settare le successive righe fino a %% -... , non modificare nient'altro!

% ATTENZIONE!! prima di usare questo script:
% 1. settare in ECLIPSE almento un parametro da identificare
% 2. assicurarsi che il modello contenga almeno una non-linearità nelle
% variabili
% 3. assicurarsi che nelle relazioni non lineari non compaiano parametri da
% identificare, se si, modificare il modello Dymola attraverso la creazione
% di variabili algebriche ausiliarie

%% ------------------------------------------------------------------------

del = [nome_m,'.mat'];
delete(del)

%--------------------------------------------------------------------------
% costruzione della funzione MATLAB 'Theta' contenente le funzioni non lineari
save('temp.mat','nome_m','nome_xml');
warning('off','all');
run(nome_m);
load temp;
delete('temp.mat');
tag = [nome_m,'_LFT'];
nomIn1 = char([nome_m, '.m']);
nomIn2 = char([nome_xml, '.xml']);
fid1 = fopen(nomIn1, 'rt');
fid2 = fopen(nomIn2, 'rt');
C1 = textscan(fid1, '%s', 'delimiter', '\n', 'whitespace', ' ', 'BufSize', 8190);
C2 = textscan(fid2, '%s', 'delimiter', '\n', 'whitespace', ' ', 'BufSize', 8190);
D1 = char(C1{:});                                                           % file .m
D2 = char(C2{:});                                                           % file .xml
clear('C1');
clear('C2');
fclose(fid1);
fclose(fid2);
clc;
        diary('command_window.mos');                                        % salvataggio temporaneo della stampa su command w.
        diary on;
        LFT
        diary off;
        id_c_w = fopen('command_window.mos', 'rt');
        C_c_w = textscan(id_c_w, '%s', 'delimiter', '\n', 'whitespace', ' ');
        D_c_w = char(C_c_w{:}); % 
        clear('C_c_w');
        fclose(id_c_w);
        delete('command_window.mos');                                       % fine salvataggio temporaneo della stampa su command w.
DimsType = {};                                                              % struttura dati contenente la struttura delle retroazioni 
nzw = 0;                                                                    % numero variabili z e w
nTypeDelta = 0;                                                             % numero di parametri Delta senza contare ripetizioni
ntheta = 0;                                                                   % numero funzioni non lineari
ind_read_struct = 1;
ind_start_size = strfind(D_c_w(strmatch('Name', D_c_w),:), 'Dims');         % indice di colonna dove inizia il titolo Dims
ind_start_type = strfind(D_c_w(strmatch('Name', D_c_w),:), 'Type');         % indice di colonna dove inizia il titolo Type
while ind_read_struct+strmatch('Name', D_c_w)<=size(D_c_w,1)                % scansione D_c_w (=copia di command w.) e costruzione della struttura dati contenente la struttura delle retroazioni
    ind_stop_size = strfind(D_c_w(strmatch('Name', D_c_w)+ind_read_struct,ind_start_size:end), ' ');
    string_size = D_c_w(strmatch('Name', D_c_w)+ind_read_struct,ind_start_size:ind_start_size+ind_stop_size(1)-2);
    type = D_c_w(strmatch('Name', D_c_w)+ind_read_struct,ind_start_type:ind_start_type+2);
    ind = strfind(D_c_w(strmatch('Name', D_c_w)+ind_read_struct,:), ' ');
    if strcmp(D_c_w(strmatch('Name', D_c_w)+ind_read_struct,ind-4:ind-1),'_par') % la stringa del primo campo è un parametro incerto
        name = D_c_w(strmatch('Name', D_c_w)+ind_read_struct,1:ind-5);
    else
        name = D_c_w(strmatch('Name', D_c_w)+ind_read_struct,1:ind-1);      % la stringa del primo campo è una funzione non lineare
    end;
    DimsType{ind_read_struct, 1} =  char(name);
    DimsType{ind_read_struct, 2} =  string_size(1:strfind(string_size, 'x')-1);
    DimsType{ind_read_struct, 3} =  string_size(strfind(string_size, 'x')+1:end);
    DimsType{ind_read_struct, 4} =  char(type);
    if strcmp(type, 'LTI')
        nzw = nzw + str2double(DimsType{ind_read_struct,2});
        nTypeDelta = nTypeDelta + 1;  
    else
        ntheta = ntheta + str2double(DimsType{ind_read_struct,2});
    end;   
    ind_read_struct = ind_read_struct +1;
end;                                                                        % fine scansione D_c_w ( = copia di command w.) e costruzione della struttura dati contenente la struttura delle retroazioni
nomOut4 = [tag, '_Theta.m'];
g4 = fopen(nomOut4, 'wt');                                                  % apro ora il file perchè stampo le funzioni non lineari durante la scansione e la sostituzione simbolica
fprintf(g4, '%s\n', ['function  y = ', tag,'_Theta(omega)']);
vett_ind_eq_nl=strmatch('% Non linear lft named', D1);                      % determina gli indici di riga -1 dove si trovano tutte le espressioni non lineari
nu = size(LFT.b,2);                                                         % numero di ingressi
nx = size(LFT.a,2)-nzw-ntheta;                                                % numero stati
nomega = size(LFT.a,1)-nx-nzw;                                              % numero ingressi funzioni non lineari
ny = size(LFT.d,1);                                                         % numero uscite
ind_u = 1;                                                                  % inizio lettura nomi delle variabili d'ingresso e di stato
ind_st = 1;
ind_ou = 1; 
elenco_ingressi = {};
elenco_stati = {};
elenco_uscite = {};
ind_riga_uxy  = strmatch('<variable id=' , D2);
for i = 1:length(ind_riga_uxy) 
    if ~isempty(strfind(D2(ind_riga_uxy(i), :),'direction="input"'))
        ind_nome_xuy = strfind(D2(ind_riga_uxy(i),:), '"');
        elenco_ingressi{ind_u} = D2(ind_riga_uxy(i),ind_nome_xuy(3)+1:ind_nome_xuy(4)-1);
        ind_u = ind_u+1;
    end;
    if ~isempty(strfind(D2(ind_riga_uxy(i), :),'variability="continuousState"'))
        ind_nome_xuy = strfind(D2(ind_riga_uxy(i),:), '"');
        elenco_stati{ind_st} = D2(ind_riga_uxy(i),ind_nome_xuy(3)+1:ind_nome_xuy(4)-1);
        ind_st = ind_st+1;
    end;
    if ~isempty(strfind(D2(ind_riga_uxy(i), :),'direction="output"'))
        ind_nome_xuy = strfind(D2(ind_riga_uxy(i),:), '"');
        elenco_uscite{ind_ou} = D2(ind_riga_uxy(i),ind_nome_xuy(3)+1:ind_nome_xuy(4)-1);
        ind_ou = ind_ou+1;
    end;
end;
vett_eq_nl_omega = {};
i = 1;
ind_omega = 1;
iii = 1;
while i <= size(DimsType,1)                                                 % scansione di ogni espressione non lineare all'interno del file.m che genera la LFT
    if strcmp(DimsType{i,4}, 'NLM')                                         % scorrendo DimsType verifica se siamo in posizione di un'eq_nl
        name_NL_eq = DimsType{i,1};                                         % legge il nome dell'equazione non lineare (Psi..)               
        ind_riga_D1_eq_nl = strmatch(['% Non linear lft named: ''',name_NL_eq,''''], D1); % cerca in D1 l'espressione con nome (Psi..)
        res = deblank(D1(ind_riga_D1_eq_nl+1,:));                           % equazione non lineare scritta per intero come residuo
        ind_meno = strfind(res, '-');                                       % prende l'indice del meno
        eq_nl = res(ind_meno(1)+1:end);                                     % estrae dal residuo l'espressione di theta
        ind_eq_nl = vett_ind_eq_nl(iii)+1;
        title=deblank(D1(ind_eq_nl-1,:));                                   
        nr_eq_nl = textscan(name_NL_eq(1,4:end), '%d');                     % numero equazione non lineare, da script per generazione LFT
        [n_omega, n_ing] = size(S{nr_eq_nl{1}+1});                          % il +1 è messo perchè la nomenclatura delle funzioni non lineari è sotto di 1
        for ii = 1:n_omega                                                  % inizio sostituzione nomi variabili
            posizione = find(S{nr_eq_nl{1}+1}(ii,:)==1);
            if posizione > nx+nu
                id_blocco_BLT = posizione-nx-nu;
                ind_blocco_xml = strmatch(['<bltBlock id="' int2str(id_blocco_BLT) '">'], D2); % mi da la riga del xml dove si associa blocco a equazione
                ind_eq_start_stop = strfind(D2(ind_blocco_xml+1,:), '"');	% tra apici trovo il numero di equazione trattata nel blocco id_blocco_BLT
                ind_eq_start = ind_eq_start_stop(1);
                ind_eq_stop = ind_eq_start_stop(2);
                id_equazione = D2(ind_blocco_xml+1,ind_eq_start+1:ind_eq_stop-1); % ma da l'id dell'equazione della BLT dove si calcola la variabile di cui mi interessa sapere il nome
                inizio_mA = strmatch('<matchingAlgorithm>', D2)+1;          % inizio descrizione matching algorithm
                fine_mA = strmatch('</matchingAlgorithm>', D2)-1;           % fine descrizione matching algorithm
                for cc = 1:fine_mA-inizio_mA+1                              % cerca la riga dove si associa id_equazione
                    ind_match_eq_var = strfind(D2(inizio_mA+cc-1,:), ['equationId="' id_equazione '" />']);
                    if ~isempty(ind_match_eq_var)                           % riga trovata
                        ind_var_start_stop = strfind(D2(inizio_mA+cc-1,:), '"'); % cerco gli indici di tutti gli apici
                        ind_var_start = ind_var_start_stop(1);              % l'id variabile è contenuto tra i primi due apici della riga
                        ind_var_stop = ind_var_start_stop(2);
                        id_variabile = D2(inizio_mA+cc-1,ind_var_start+1:ind_var_stop-1);
                    end;
                end;
                ind_var_xml = strmatch(['<variable id="' id_variabile '"'],D2); % mi da due indici, il primo è una var alg o di stato o di alias, il secondo è una u
                ind_nome_start_stop = strfind(D2(ind_var_xml(1),:), '"');   % prendo il primo indice perchè ho bisogno di una var alg
                ind_nome_start = ind_nome_start_stop(3)+1;                  % indice inizio nome variabile algebrica, si trova tra i secondi 2 apici
                ind_nome_stop = ind_nome_start_stop(4)-1;                   % indice fine nome
                nome_var = D2(ind_var_xml(1),ind_nome_start:ind_nome_stop);
            else
                if posizione > nx %% è una u
                    id_input = posizione-nx;
                    nome_var = elenco_ingressi{id_input};           
                else   %% è uno stato
                    nome_var = elenco_stati{posizione};
                end;
            end;
            eq_nl = regexprep(eq_nl, nome_var, ['omega(' int2str(ind_omega) ')']); % sostituzione simbolica
            ind_omega = ind_omega+1;
        end;                                                                % fine sostituzione nomi variabili
        indper = strfind(res(1:ind_meno-1), '*');
        if ~isempty(indper)                                                 % controlla se prima del meno la theta è stata moltiplicata per qualche fattore
            vett_eq_nl_omega{iii,1} = char([eq_nl, '/',res(3:indper-1)]);     % copia l'espressione della funzione non lineare divisa per il fatt moltiplicante
        else
            vett_eq_nl_omega{iii,1} = eq_nl;                                  % copia l'espressione della funzione non lineare così come l'ha letta finora
        end;
        vett_eq_nl_omega{iii,2} = n_omega;
        iii=iii+1;
    end;
    i = i+1; 
end;
ind_omega = ind_omega-1;
par_nl = {};                                                                % struttura dati per i parametri delle funzioni non lineari
pp = 1;
indd = strmatch('<variable id="', D2);
for cc = 1:length(indd)                                                     % inizio ricerca parametri generici nell'xml
    if strfind(D2(indd(cc),:), 'parameter')                                 % trovato una riga contenente un parametro
        in_fin = strfind(D2(indd(cc),:), '"');
        nome_parametro = char(D2(indd(cc),in_fin(3)+1:in_fin(4)-1));        % leggo il nome del parametro
        kk = 1;
        trovato = false;
        while (~trovato)&&(kk<=size(vett_eq_nl_omega,1))                    % inizio ricerca del parametro 'nome_parametro' nelle funzioni non lineari
            if ~isempty(strfind(vett_eq_nl_omega{kk,1}, nome_parametro))&&(trovato==false);  % il parametro coincide con un parametro delle funzioni non lienari
                par_nl{pp} = nome_parametro;                                % riempimento della struttura dati per i parametri delle funzioni non lineari
                pp = pp+1;
                trovato = true;                                             % il parametro è un parametro delle funzioni non lineari
            end;
            kk = kk+1;
        end;                                                                % ricerca del parametro 'nome_parametro' nelle funzioni non lineari
    end;
end;                                                                       % fine ricerca parametri generici

NLEquationsParValues = struct();
ind_riga = strmatch('<variable id="',D2);
for i=1:length(par_nl)
    for ii = 1:length(ind_riga)
        if ~isempty(strfind(D2(ind_riga(ii), :), ['"',par_nl{i},'"']))
            if ~isempty(strfind(D2(ind_riga(ii)+2, :), 'string'))
                ind_in_fin_valore = strfind(D2(ind_riga(ii,:)+2,:), '"');
                NLEquationsParValues.(par_nl{i}) = D2(ind_riga(ii,:)+2,ind_in_fin_valore(1)+1:ind_in_fin_valore(2)-1);
            else
                ind_in_fin_valore = strfind(D2(ind_riga(ii,:)+3,:), '"');
                NLEquationsParValues.(par_nl{i}) = D2(ind_riga(ii,:)+3,ind_in_fin_valore(1)+1:ind_in_fin_valore(2)-1);
            end;            
        end;
    end;
end;

    %--------------------------------------------------------------------------
    % costruzione della funzione MATLAB 'dThetadOmega' contenente le derivate delle funzioni non lineari rispetto ai loro argomenti 
        % costruzione temporanea script per derivazione simbolica
        nomOut8 = 'temp.m';
        g8 = fopen(nomOut8, 'wt');
        for i = 1:ind_omega
            fprintf(g8, '%s\n', ['syms omega_',int2str(i)]);
        end;
        i = 1;
        par_nl_underscore = par_nl;
        while i <= length(par_nl)
            par_nl_underscore{i} = regexprep(par_nl{i}, '\.', '_'); % sostituzione simbolica
            fprintf(g8, '%s\n', ['syms ', par_nl_underscore{i}]);
            i=i+1; 
        end;
        fprintf(g8, '\n');   
        fprintf(g8, '%s\n', 'omega = [');
        for i = 1:ind_omega
            fprintf(g8, '%s\n', ['omega_',int2str(i),';']);
        end;
        fprintf(g8, '%s\n', '];'); 
        fprintf(g8, '%s\n', 'y =[');
        vett_eq_nl_omega_underscore = vett_eq_nl_omega;
        for i = 1:length(vett_ind_eq_nl)
            vett_eq_nl_omega_underscore{i,1} = vett_eq_nl_omega{i,1};
            for j = 1:length(par_nl)
                vett_eq_nl_omega_underscore{i,1} = regexprep(vett_eq_nl_omega_underscore{i,1}, par_nl{j}, par_nl_underscore{j});
            end;    
            fprintf(g8, '    %s\n', vett_eq_nl_omega_underscore{i,1});  
        end;
        fprintf(g8, '%s\n', '];');

        fclose(g8);
        run('temp')
        testo_jacobiano = jacobian(y, omega);
        delete('temp.m');
        % costruzione temporanea script per derivazione simbolica
        nomOut7 = [tag, '_dThetadOmega.m'];
        g7 = fopen(nomOut7, 'wt');
        fprintf(g7, '%s\n', ['function y = ', tag,'_dThetadOmega(omega)']);
        i = 1;
        while i <= length(par_nl)
            fprintf(g4, '%s\n', ['  ', char(par_nl_underscore{i}), ' = ', NLEquationsParValues.(char(par_nl_underscore{i})), ';']);
            fprintf(g7, '%s\n', ['  ', char(par_nl_underscore{i}), ' = ', NLEquationsParValues.(char(par_nl_underscore{i})), ';']);
            i=i+1; 
        end;
        for i = 1:ind_omega
            fprintf(g7, '%s\n', ['omega_',int2str(i),'=omega(',int2str(i),');']);
        end;
        fprintf(g7, '%s\n', 'y=[');
        for i = 1:size(testo_jacobiano,1)
            if i == 1 
                for j = 1:size(testo_jacobiano,2)-1
                    fprintf(g7, '%s', char(testo_jacobiano(i,j)), ',    ');       
                end; %fine colonne
                fprintf(g7, '%s\n', [char(testo_jacobiano(i,size(testo_jacobiano,2))), ';']);
            else
                for j = 1:size(testo_jacobiano,2)-1
                    if j ==1
                        fprintf(g7, '%s', ['      ',char(testo_jacobiano(i,j)), ',    ']);
                    else
                        fprintf(g7, '%s', char(testo_jacobiano(i,j)), ',    ');
                    end;
                end; %fine colonne
                fprintf(g7, '%s\n', [char(testo_jacobiano(i,size(testo_jacobiano,2))), ';']);
            end;
        end; % fine righe
        fprintf(g7, '%s\n', '];');
        fclose(g7);
fprintf(g4, '%s\n', 'y =[');
for i = 1:length(vett_ind_eq_nl)
     fprintf(g4, '    %s\n', vett_eq_nl_omega_underscore{i,1});  
end;
fprintf(g4, '%s\n', '];');
fclose(g4);


%--------------------------------------------------------------------------
% riordino della forma LTI
A = LFT.a(1:nx,1:nx);
B1 = zeros(nx,nzw);
B2 = zeros(nx,ntheta);
B3 = LFT.b(1:nx,:);
C1 = zeros(nzw,nx);
D11 = zeros(nzw,nzw);
D12 = zeros(nzw,ntheta);
D13 = zeros(nzw,nu);
C2 = zeros(nomega,nx);
D21 = zeros(nomega,nzw);
D22 = zeros(nomega,ntheta);
D23 = zeros(nomega,nu);
C3 = LFT.c(:,1:nx);
D31 = zeros(size(LFT.c,1),nzw);
D32 = zeros(size(LFT.c,1),ntheta);
D33 = LFT.d;
ind_r_D11 = 1;
ind_r_D12 = 1;
ind_r_D21 = 1;
ind_r_D22 = 1;
n = ind_read_struct-1;
m = ind_read_struct-1;
puntatore_r = nx+1;
for ii=1:n
    puntatore_c = nx+1;
    if strcmp(DimsType{ii,4}, 'LTI')                                        % in riga ho una z allora riempio D11 o D12
        ind_c_D11 = 1;
        ind_c_D12 = 1;
        for jj=1:m
            if strcmp(DimsType{jj,4}, 'LTI')                                % in colonna ho una z allora riempio D11                
                D11(ind_r_D11:ind_r_D11+str2double(DimsType{ii,3})-1,ind_c_D11:ind_c_D11+str2double(DimsType{jj,2})-1)=...
                    LFT.a(puntatore_r:puntatore_r+str2double(DimsType{ii,3})-1,puntatore_c:puntatore_c+str2double(DimsType{jj,2})-1);
                ind_c_D11=ind_c_D11+str2double(DimsType{jj,2});                   
            else                                                            % in colonna ho una theta allora riempio D12           
                D12(ind_r_D12:ind_r_D12+str2double(DimsType{ii,3})-1,ind_c_D12:ind_c_D12+str2double(DimsType{jj,2})-1)=...
                    LFT.a(puntatore_r:puntatore_r+str2double(DimsType{ii,3})-1,puntatore_c:puntatore_c+str2double(DimsType{jj,2})-1);
                ind_c_D12=ind_c_D12+str2double(DimsType{jj,2});   
            end;
            puntatore_c=puntatore_c+str2double(DimsType{jj,2});
        end;        
        C1(ind_r_D11:ind_r_D11+str2double(DimsType{ii,3})-1,:) =...
            LFT.a(puntatore_r:puntatore_r+str2double(DimsType{ii,3})-1,1:nx);
        D13(ind_r_D11:ind_r_D11+str2double(DimsType{ii,3})-1,:) =...
            LFT.b(puntatore_r:puntatore_r+str2double(DimsType{ii,3})-1,:);
        ind_r_D11=ind_r_D11+str2double(DimsType{ii,3});
        ind_r_D12=ind_r_D12+str2double(DimsType{ii,3});
    else                                                                    % in riga ho una omega allora riempio D21 o D22
        ind_c_D21 = 1;
        ind_c_D22 = 1;
        for jj=1:m
            if strcmp(DimsType{jj,4}, 'LTI')                                % in colonna ho una z allora riempio D21                
                D21(ind_r_D21:ind_r_D21+str2double(DimsType{ii,3})-1,ind_c_D21:ind_c_D21+str2double(DimsType{jj,2})-1)=...
                LFT.a(puntatore_r:puntatore_r+str2double(DimsType{ii,3})-1,puntatore_c:puntatore_c+str2double(DimsType{jj,2})-1);
                ind_c_D21=ind_c_D21+str2double(DimsType{jj,2});                   
            else                                                            % in colonna ho una theta allora riempio D22           
                D22(ind_r_D22:ind_r_D22+str2double(DimsType{ii,3})-1,ind_c_D22:ind_c_D22+str2double(DimsType{jj,2})-1)=...
                LFT.a(puntatore_r:puntatore_r+str2double(DimsType{ii,3})-1,puntatore_c:puntatore_c+str2double(DimsType{jj,2})-1);
                ind_c_D22=ind_c_D22+str2double(DimsType{jj,2});   
            end;
            puntatore_c=puntatore_c+str2double(DimsType{jj,2});
        end;
        C2(ind_r_D21:ind_r_D21+str2double(DimsType{ii,3})-1,:) =...
            LFT.a(puntatore_r:puntatore_r+str2double(DimsType{ii,3})-1,1:nx);
        D23(ind_r_D21:ind_r_D21+str2double(DimsType{ii,3})-1,:) =...
            LFT.b(puntatore_r:puntatore_r+str2double(DimsType{ii,3})-1,:);
        ind_r_D21=ind_r_D21+str2double(DimsType{ii,3});
        ind_r_D22=ind_r_D22+str2double(DimsType{ii,3});
    end;
    puntatore_r=puntatore_r+str2double(DimsType{ii,3});
end;
puntatore_c = nx+1;
ind_c_B1 = 1;
ind_c_B2 = 1;
for jj=1:m
    if strcmp(DimsType{jj,4}, 'LTI')
        B1(:,ind_c_B1:ind_c_B1+str2double(DimsType{jj,2})-1)=...
            LFT.a(1:nx,puntatore_c:puntatore_c+str2double(DimsType{jj,2})-1);
        D31(:,ind_c_B1:ind_c_B1+str2double(DimsType{jj,2})-1)=...
             LFT.c(:,puntatore_c:puntatore_c+str2double(DimsType{jj,2})-1);
        ind_c_B1=ind_c_B1+str2double(DimsType{jj,2});                   
    else
        B2(:,ind_c_B2:ind_c_B2+str2double(DimsType{jj,2})-1)=...
            LFT.a(1:nx,puntatore_c:puntatore_c+str2double(DimsType{jj,2})-1);
        D32(:,ind_c_B2:ind_c_B2+str2double(DimsType{jj,2})-1)=...
             LFT.c(:,puntatore_c:puntatore_c+str2double(DimsType{jj,2})-1);
        ind_c_B2=ind_c_B2+str2double(DimsType{jj,2});
    end;
    puntatore_c=puntatore_c+str2double(DimsType{jj,2});
end;   
LTI = struct('A',A,     'B1',B1,    'B2',B2,    'B3',B3,...
             'C1',C1,   'D11',D11,  'D12',D12,  'D13',D13,...
             'C2',C2,   'D21',D21,  'D22',D22,  'D23',D23,...
             'C3',C3,   'D31',D31,  'D32',D32,  'D33',D33);


%--------------------------------------------------------------------------
% info blocchi retroazione  
ind_DimsType=1;
zw=1;
ind_nome_par_inc=1;
v_DELTA = sym(zeros(nzw,1));
delta_min = {};
delta_max = {};
list_TypeDelta = {};
ind_list_TypeDelta = 1;
ind_min_max = 1;
ind_quadre_chiuse = strfind(D1(strmatch('lfrs', D1),:),']');
ind_spazio_par = strfind(D1(strmatch('lfrs', D1),1:ind_quadre_chiuse(2)+1), ' ');

while ind_DimsType<=size(DimsType,1)
    if strcmp(DimsType{ind_DimsType,4}, 'LTI')
        list_TypeDelta{ind_list_TypeDelta, 1} = DimsType{ind_DimsType,1};
        list_TypeDelta{ind_list_TypeDelta, 2} = DimsType{ind_DimsType,2};
        ind_list_TypeDelta = ind_list_TypeDelta+1;
        n_ricorrenze = str2double(DimsType{ind_DimsType,2});
        for rr=1:n_ricorrenze
            v_DELTA(zw)=DimsType{ind_DimsType,1};
            zw=zw+1;
        end;
        ind_parMATLAB = strfind(D1(strmatch('lfrs', D1),:),DimsType{ind_DimsType,1})-1;
        ord_parMATLAB = find(ind_spazio_par == ind_parMATLAB);
        
        if D1(strmatch('lfrs', D1),ind_spazio_par(nTypeDelta+ord_parMATLAB)+1)=='['
            if D1(strmatch('lfrs', D1),ind_spazio_par(nTypeDelta+ord_parMATLAB+1)-1)==']'
                delta_min{ind_min_max} = D1(strmatch('lfrs', D1),ind_spazio_par(nTypeDelta+ord_parMATLAB)+2:ind_spazio_par(nTypeDelta+ord_parMATLAB+1)-2);
                delta_max{ind_min_max} = D1(strmatch('lfrs', D1),ind_spazio_par(2*nTypeDelta+ord_parMATLAB)+2:ind_spazio_par(2*nTypeDelta+ord_parMATLAB+1)-2);
            else
                delta_min{ind_min_max} = D1(strmatch('lfrs', D1),ind_spazio_par(nTypeDelta+ord_parMATLAB)+2:ind_spazio_par(nTypeDelta+ord_parMATLAB+1)-1);
                delta_max{ind_min_max} = D1(strmatch('lfrs', D1),ind_spazio_par(2*nTypeDelta+ord_parMATLAB)+2:ind_spazio_par(2*nTypeDelta+ord_parMATLAB+1)-1);
            end;
        else
            if D1(strmatch('lfrs', D1),ind_spazio_par(nTypeDelta+ord_parMATLAB+1)-1)==']'
                delta_min{ind_min_max} = D1(strmatch('lfrs', D1),ind_spazio_par(nTypeDelta+ord_parMATLAB)+1:ind_spazio_par(nTypeDelta+ord_parMATLAB+1)-2);
                delta_max{ind_min_max} = D1(strmatch('lfrs', D1),ind_spazio_par(2*nTypeDelta+ord_parMATLAB)+1:ind_spazio_par(2*nTypeDelta+ord_parMATLAB+1)-2);
            else
                delta_min{ind_min_max} = D1(strmatch('lfrs', D1),ind_spazio_par(nTypeDelta+ord_parMATLAB)+1:ind_spazio_par(nTypeDelta+ord_parMATLAB+1)-1);
                delta_max{ind_min_max} = D1(strmatch('lfrs', D1),ind_spazio_par(2*nTypeDelta+ord_parMATLAB)+1:ind_spazio_par(2*nTypeDelta+ord_parMATLAB+1)-1);
            end;
        end;
        ind_min_max = ind_min_max+1;
    end;
    ind_DimsType=ind_DimsType+1;
end;
DELTA = sym(diag(v_DELTA));


%--------------------------------------------------------------------------
% analisi di sparsità
A_Sparsity = sparsityMatrix(A);
B1_Sparsity = sparsityMatrix(B1);
B2_Sparsity = sparsityMatrix(B2);
B3_Sparsity = sparsityMatrix(B3);
C1_Sparsity = sparsityMatrix(C1);
C2_Sparsity = sparsityMatrix(C2);
C3_Sparsity = sparsityMatrix(C3);
D11_Sparsity = sparsityMatrix(D11);
D12_Sparsity = sparsityMatrix(D12);
D13_Sparsity = sparsityMatrix(D13);
D21_Sparsity = sparsityMatrix(D21);
D22_Sparsity = sparsityMatrix(D22);
D23_Sparsity = sparsityMatrix(D23);
D31_Sparsity = sparsityMatrix(D31);
D32_Sparsity = sparsityMatrix(D32);
D33_Sparsity = sparsityMatrix(D33);
dThetadOmega_Sparsity = true(size(testo_jacobiano));

[row_zero_elements, col_zero_elements] = find(testo_jacobiano==0);
if ~isempty(find(testo_jacobiano==0))
    for ind_coord = 1:size(row_zero_elements,1)
        dThetadOmega_Sparsity(row_zero_elements(ind_coord),col_zero_elements(ind_coord))=false;
    end;
end;

VCT8_Sparsity = rowColProductBoolean(invBoolean([((diag(true(size(C1,1),1)))|D11_Sparsity), rowColProductBoolean(D12_Sparsity,dThetadOmega_Sparsity)
                                                 D21_Sparsity,                              ((diag(true(size(C2,1),1)))|rowColProductBoolean(D22_Sparsity,dThetadOmega_Sparsity))]),[C1_Sparsity, D13_Sparsity, D11_Sparsity
                                                                                                                                                                                     C2_Sparsity, D23_Sparsity, D21_Sparsity]);                                                 
VCT9_Sparsity = ([A_Sparsity,  B1_Sparsity
                  C3_Sparsity, D31_Sparsity]|rowColProductBoolean([B1_Sparsity,  rowColProductBoolean(B2_Sparsity,dThetadOmega_Sparsity)
                                                                   D31_Sparsity, rowColProductBoolean(D32_Sparsity,dThetadOmega_Sparsity)],[VCT8_Sparsity(:,1:size(A,2)),VCT8_Sparsity(:,size(A,2)+size(D13,2)+1:end)]));

VCT10_Sparsity = ([A_Sparsity, B3_Sparsity
                  C3_Sparsity, D33_Sparsity]|rowColProductBoolean([B1_Sparsity,  rowColProductBoolean(B2_Sparsity,dThetadOmega_Sparsity)
                                                                   D31_Sparsity, rowColProductBoolean(D32_Sparsity,dThetadOmega_Sparsity)],VCT8_Sparsity(:,1:size(A,2)+size(D13,2))));
                                      

SparsityInformations = struct();


[row,col] = find(VCT9_Sparsity(1:size(A,1),1:size(A,2)));
SparsityInformations.VCT9{1,1}.NumNonZeroElements = length(row);
SparsityInformations.VCT9{1,1}.row = row;
SparsityInformations.VCT9{1,1}.col = col;

[row,col] = find(VCT9_Sparsity(size(A,1)+1:end,1:size(A,2)));
SparsityInformations.VCT9{2,1}.NumNonZeroElements = length(row);
SparsityInformations.VCT9{2,1}.row = row;
SparsityInformations.VCT9{2,1}.col = col;

ind_stop_z = 0;
for i=1:size(list_TypeDelta,1)
    ind_start_z = ind_stop_z+1;
    ind_stop_z = ind_start_z+str2double(list_TypeDelta{i,2})-1;
    list_TypeDelta{i,3} = ind_start_z;
    list_TypeDelta{i,4} = ind_stop_z;

    [row,col] = find(VCT9_Sparsity(1:size(A,1),size(A,2)+ind_start_z:size(A,2)+ind_stop_z));
    SparsityInformations.VCT9{1,1+i}.NumNonZeroElements = length(row);
    SparsityInformations.VCT9{1,1+i}.row = row;
    SparsityInformations.VCT9{1,1+i}.col = col;
    
    [row,col] = find(VCT9_Sparsity(size(A,1)+1:end,size(A,2)+ind_start_z:size(A,2)+ind_stop_z));
    SparsityInformations.VCT9{2,1+i}.NumNonZeroElements = length(row);
    SparsityInformations.VCT9{2,1+i}.row = row;
    SparsityInformations.VCT9{2,1+i}.col = col;
end;

[row,col] = find(VCT10_Sparsity(1:size(A,1),1:size(A,2)));
SparsityInformations.VCT10{1,1}.NumNonZeroElements = length(row);
SparsityInformations.VCT10{1,1}.row = row;
SparsityInformations.VCT10{1,1}.col = col;

[row,col] = find(VCT10_Sparsity(size(A,1)+1:end,1:size(A,2)));
SparsityInformations.VCT10{2,1}.NumNonZeroElements = length(row);
SparsityInformations.VCT10{2,1}.row = row;
SparsityInformations.VCT10{2,1}.col = col;


for i=1:size(B3,2)
    [row,col] = find(VCT10_Sparsity(1:size(A,1),size(A,2)+i));
    SparsityInformations.VCT10{1,1+i}.NumNonZeroElements = length(row);
    SparsityInformations.VCT10{1,1+i}.row = row;
    SparsityInformations.VCT10{1,1+i}.col = col;

    [row,col] = find(VCT10_Sparsity(size(A,1)+1:end,size(A,2)+i));
    SparsityInformations.VCT10{2,1+i}.NumNonZeroElements = length(row);
    SparsityInformations.VCT10{2,1+i}.row = row;
    SparsityInformations.VCT10{2,1+i}.col = col;
end;



DeltaSym = {'parName', 'indStartDiag', 'indStopDiag', 'LowerBound', 'HigherBound'};
for pp=1:size(list_TypeDelta,1)
    DeltaSym = cat(1, DeltaSym,{list_TypeDelta{pp,1},list_TypeDelta{pp,3},list_TypeDelta{pp,4}, str2double(delta_min{pp}), str2double(delta_max{pp})} );
end;

% DimsType
lftfun = struct('LTI', LTI,...
                'DeltaSym', {DeltaSym},...
                'DeltaVal', zeros(nzw),...
                'Theta', eval(['@', tag, '_Theta']),...
                'dThetadOmega', eval(['@', tag, '_dThetadOmega']),...
                'Sparsity', SparsityInformations,...
                'ISO', struct('InputVariablesNames', {elenco_ingressi},...
                              'StateVariablesNames', {elenco_stati},...
                              'OutputVariablesNames', {elenco_uscite},...
                              'Nominal', ones(1,size(C3,1))));               
                  
eval([tag, ' = lftfun']);
save([tag, '.mat'], tag); 

clear all;
clc;
warning('on','all');

'done!'










