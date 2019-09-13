function [ D ] = decimation( M, dec )
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

[r c]=size(M);
D=[];
for i=1:fix(r/dec)
    D(i,:)=M(i*dec,:);
end
end

