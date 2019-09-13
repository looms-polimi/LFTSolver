function [index,direction]=minFind(M)
index_r=1;  % index_r alla fine indicherà la riga in cui ci sono più zeri
nr=length(M(:,1));
lmin_row=nr;
for i=1:nr  % scansiono tutte le righe
    temp=find(M(i,:));
    if isempty(temp)    % se è vuota è quella con più zeri
        lmin_row=0;
        index_r=i;
        break;
    end;
    if length(temp)<lmin_row  % se il numero di posti non zero è minore del corrente
        lmin_row=length(temp);
        index_r=i;
    end;
end;
index_c=1;  % index_c alla fine indicherà la riga in cui ci sono più zeri
nc=length(M(1,:));
lmin_col=nc;
for i=1:nc  % scansiono tutte le colonne
    temp=find(M(:,i));
    if isempty(temp)  % se è vuota è quella con più zeri
        lmin_col=0;
        index_c=i;
        break;
    end;
    if length(temp)<lmin_col  % se il numero di posti non zero è minore del corrente
        lmin_col=length(temp);
        index_c=i;
    end;
end;
lmin=min(lmin_row,lmin_col);
if lmin==lmin_row
    direction='row';
    index=index_r;
else
    direction='col';
    index=index_c;
end;
end
