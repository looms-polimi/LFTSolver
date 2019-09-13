function M_red=reduceMat(M,row,col)
[nr,nc]=size(M);
if row==1 && col==1
    M_red=M(row+1:nr,col+1:nc);
    return;
end;
if row==1  % ma col è maggiore di 1
    M_left=M(row+1:nr,1:col-1);
    if col<nc
        M_right=M(row+1:nr,col+1:nc);
    else % col==nc
        M_right=[];
    end;
    M_red=horzcat(M_left,M_right);
    return;
end;
if col==1 % ma row è maggiore di 1
    M_up=M(1:row-1,col+1:nc);
    if row<nr
        M_down=M(row+1:nr,col+1:nc);
    else % row=nr
        M_down=[];
    end;
    M_red=vertcat(M_up,M_down);
    return;
end;
% row è maggiore di 1 e col è maggiore di 1
M_up_left=M(1:row-1,1:col-1);
if col<nc
    M_up_right=M(1:row-1,col+1:nc);
else 
    M_up_right=[];
end;
M_up=horzcat(M_up_left,M_up_right);
if row<nr
    M_down_left=M(row+1:nr,1:col-1);
    if col<nc
        M_down_right=M(row+1:nr,col+1:nc);
    else 
        M_down_right=[];
    end;
else
    M_down_left=[];
    M_down_right=[];
end;
M_down=horzcat(M_down_left,M_down_right);
M_red=vertcat(M_up,M_down);
end