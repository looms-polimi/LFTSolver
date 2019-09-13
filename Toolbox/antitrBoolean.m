function res=antitrBoolean(M)
% if not(islogical(M))
%     error ('Only logical matrices are accepted. Use det instead');
% end;
if length(M(:,1))~=length(M(1,:))
    error ('Matrix must be square');
end;
res=true;
ad=fliplr(diag(M));
I=find(ad);
if length(I)~=length(ad)
    res=false;
end;
end