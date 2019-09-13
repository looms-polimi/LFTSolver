function res=trBoolean(M)
% if not(islogical(M))
%     error ('Only logical matrices are accepted. Use tr instead');
% end;
if length(M(:,1))~=length(M(1,:))
    error ('Matrix must be square');
end;
res=true;
d=diag(M);
I=find(d);
if length(I)~=length(d)
    res=false;
end;
end