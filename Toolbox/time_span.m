function [ vo ] = time_span( v )
% Allargo il time input a ore

temp=zeros(length(v)*2-1,1);
imax=length(v)-1;
for i=1:imax
    temp(2*i-1,1)=v(i,1);
    temp(2*i)=(v(i,1)+v(i+1,1))/2;
end
temp(length(v)*2-1,1)=v(end,1);
vo=temp;

