% Product Row Column Booleans
function y=rowColProductBoolean(M1,M2)
% Computes the product Row Columns between Boolean matrices
if not(islogical(M1) && islogical(M2))
    M1=logical(M1);
    M2=logical(M2);
end;
[m,n1]=size(M1);
[n2,p]=size(M2);
if not(isequal(n1,n2))
    error ('Non comformable matrices');
end;
n=n1; 
if isequal([m,n],[1,1]) && isequal([n,p],[1,1])
    y=M1&&M2;
    return;
end;
y=false(m,p);
for i=1:m  % per ogni riga della prima matrice 
    if isequal(M1(i,:),false(1,n))
        y(i,:)=false(1,p);
        continue;  % salta alla prossima riga
    end;
    for j=1:p  % per ogni colonna della seconda
        if isequal(M2(:,j),false(n,1))
            y(:,j)=false(m,1);
            continue;  % salta alla prossima colonna
        end;
        if y(i,j)==true
            continue;
        end;
        for k=1:n 
            if y(i,j)==true % appena è TRUE salta al successivo
                continue;
            end;
            y(i,j)=y(i,j)||(M1(i,k)&& M2(k,j));
        end;
    end;
end;
end