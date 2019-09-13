function y=detBoolean(M)
% Calcoliamo il determinante usando la regola seguente:
%   ogni moltiplicazione diventa un &&
%   ogni somma diventa un ||

% if not(islogical(M))
%     error ('Only logical matrices are accepted. Use det instead');
% end;
if length(M(:,1))~=length(M(1,:))
    error ('Matrix must be square');
end;
%% casi base:
if length(M(1,:))==1     % dim=1
    y=M;
    return;
end;
if isequal(triu(M),M)||isequal(tril(M),M)  % triangular 
    y=trBoolean(M);
    return;
end;
if trBoolean(M)==true
    y=true;
    return;
end;
if antitrBoolean(M)==true
    y=true;
    return;
end;
if length(M(1,:))==2     % dim=2
    y=(M(1,1)&&M(2,2))||(M(1,2)&&M(2,1));
    return;
end;
if length(M(1,:))==3     % dim=3
    [ind,line]=minFind(M); % trova la riga/colonna dove ci sono più zeri
    if isequal(line,'row')
        inds=find(M(ind,:));
    else
        inds=find(M(:,ind));
    end;
    if isempty(inds)
        y=false;
        return;
    end;
    y=false;
    for i=1:length(inds)
        if y==true
            return;
        end;
        if isequal(line,'row')
            M_red=reduceMat(M,ind,inds(i));
        else
            M_red=reduceMat(M,inds(i),ind);
        end;
        y=y||((M_red(1,1)&&M_red(2,2))||((M_red(1,2)&&M_red(2,1))));
    end;
    return;
end;
%% passo induttivo (n>=4)
[ind,line]=minFind(M);
if isequal(line,'row')
    inds=find(M(ind,:));
else
    inds=find(M(:,ind));
end;
if isempty(inds)
    y=false;
    return;
end;
y=false;
for i=1:length(inds)
    if y==true
        return;
    end;
    if isequal(line,'row')
        M_red=reduceMat(M,ind,inds(i));
    else
        M_red=reduceMat(M,inds(i),ind);
    end;
    y=y||detBoolean(M_red);  %chiamata ricorsiva
end;

end
