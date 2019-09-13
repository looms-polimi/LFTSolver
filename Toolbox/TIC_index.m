function [ val ] = TIC_index( y, ym )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if size(y)==size(ym)
    display('Computing Theil’s inequality coefficient (TIC; Theil, 1961)');
else
    display('Error, invalid input dimension');
    exit;
end

sum1=0;
sum2=0;
sum3=0;
for i=1:length(y)
    if isnan(y(i)) || isnan(ym(i))
        
    else
        sum1=sum1+(y(i)-ym(i))^2;
        sum2=sum2+y(i)^2;
        sum3=sum3+ym(i)^2;
    end
end

val=sqrt(sum1)/(sqrt(sum2)+sqrt(sum3));
%display ('TIC = ', val);

end

