function Minv=invBoolean(M)
if detBoolean(M)==false
    error('The matrix is not invertible!');
end;
[nr,nc]=size(M);
if not(nr==nc)
    error ('The matrix must be square!');
end;
n=nr;
if n==1
    Minv=M;
    return;
end;
Minv=false(n);
for i=1:n
    for j=1:n
        Mred=reduceMat(M,i,j);
        Minv(i,j)=detBoolean(Mred);
    end;
end;
Minv=Minv';
end