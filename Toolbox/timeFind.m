function [y1,y2] = timeFind(t,t0,k)
    n = length(t);
    if n==2
        y1 = k+1;
        y2 = k+2;
    else
        i = floor(n/2);
        if t0>=t(i)
            if t0<t(i+1) %=
                y1 = k+i;
                y2 = k+i+1;
            else
                [y1,y2] = timeFind(t(i+1:end),t0,k+i);
            end;
        else  
            if t0>=t(i-1)
                y1 = k+i-1;
                y2 = k+i;
            else
                [y1,y2] = timeFind(t(1:i-1),t0,k);
            end;
        end;
    end;
end