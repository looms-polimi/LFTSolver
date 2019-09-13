function [lftfun] = lft_finalize(lftfun)
%% create theta and dtheta from symbolic formulas
% theta
text = '@(om) [';
for i = 1 : length(lftfun.theta_sym)
   f = char(lftfun.theta_sym(i));
   for n = 1 : lftfun.om_count
       f = strrep(f, [ 'om' num2str(n) '_1' ], [ 'om(' num2str(n) ')' ]);
   end
   text = [ text f ';' ];
end
lftfun.Theta = eval([ text ']' ]);

% dtheta
text = '@(om) [';
for i = 1 : length(lftfun.theta_sym)
   for j = 1 : lftfun.om_count
       f = char(diff(lftfun.theta_sym(i), [ 'om' num2str(j) '_1' ]));
       for n = 1 : lftfun.om_count
           f = strrep(f, [ 'om' num2str(n) '_1' ], [ 'om(' num2str(n) ')' ]);
       end
       text = [ text f ];
       if j ~= lftfun.om_count
           text = [ text ',' ];
       end
   end
   text = [ text ';' ];
end
lftfun.dThetadOmega = eval([ text ']' ]);

%% normalize deltas
E = [];
F = [];
dim=size(lftfun.DeltaSym);
jmax=dim(1,1)-1;
for j = 1 : jmax
    imin=cell2mat(lftfun.DeltaSym(1+j,2));
    imax=cell2mat(lftfun.DeltaSym(1+j,3));
    for i = imin : imax
        dmin = cell2mat(lftfun.DeltaSym(1+j,4));
        dmax = cell2mat(lftfun.DeltaSym(1+j,5));
        E = [E (dmin + dmax)/2];
        F = [F (dmax - dmin)/2];
    end
end

E = diag(E);
F = diag(F);

A=lftfun.LTI.A;
B1=lftfun.LTI.B1;
B2=lftfun.LTI.B2;
B3=lftfun.LTI.B3;
C1=lftfun.LTI.C1;
C2=lftfun.LTI.C2;
C3=lftfun.LTI.C3;
D11=lftfun.LTI.D11;
D12=lftfun.LTI.D12;
D13=lftfun.LTI.D13;
D21=lftfun.LTI.D21;
D22=lftfun.LTI.D22;
D23=lftfun.LTI.D23;
D31=lftfun.LTI.D31;
D32=lftfun.LTI.D32;
D33=lftfun.LTI.D33;

G = inv(eye(size(E))-D11*E);

lftfun.LTI.A=A+B1*E*G*C1;
lftfun.LTI.B1=B1+B1*E*G*D11;
lftfun.LTI.B2=B2+B1*E*G*D12;
lftfun.LTI.B3=B3+B1*E*G*D13;
lftfun.LTI.C1=F*G*C1;
lftfun.LTI.C2=C2+D21*E*G*C1;
lftfun.LTI.C3=C3+D31*E*G*C1;
lftfun.LTI.D11=F*G*D11;
lftfun.LTI.D12=F*G*D12;
lftfun.LTI.D13=F*G*D13;
lftfun.LTI.D21=D21+D21*E*G*D11;
lftfun.LTI.D22=D22+D21*E*G*D12;
lftfun.LTI.D23=D23+D21*E*G*D13;
lftfun.LTI.D31=D31+D31*E*G*D11;
lftfun.LTI.D32=D32+D31*E*G*D12;
lftfun.LTI.D33=D33+D31*E*G*D13;

%% generate sparsity matrices
% based on lftMatlabModel

oms = [];
for n = 1 : size(lftfun.LTI.C2,1)
    oms = [ oms; sym([ 'om' num2str(n) '_1' ]) ];
end
testo_jacobiano = jacobian(lftfun.theta_sym, oms);

if all(testo_jacobiano == 0)
    error('jacobiano nullo')
end

% analisi di sparsità
A_Sparsity = sparsityMatrix(lftfun.LTI.A);
B1_Sparsity = sparsityMatrix(lftfun.LTI.B1);
B2_Sparsity = sparsityMatrix(lftfun.LTI.B2);
B3_Sparsity = sparsityMatrix(lftfun.LTI.B3);
C1_Sparsity = sparsityMatrix(lftfun.LTI.C1);
C2_Sparsity = sparsityMatrix(lftfun.LTI.C2);
C3_Sparsity = sparsityMatrix(lftfun.LTI.C3);
D11_Sparsity = sparsityMatrix(lftfun.LTI.D11);
D12_Sparsity = sparsityMatrix(lftfun.LTI.D12);
D13_Sparsity = sparsityMatrix(lftfun.LTI.D13);
D21_Sparsity = sparsityMatrix(lftfun.LTI.D21);
D22_Sparsity = sparsityMatrix(lftfun.LTI.D22);
D23_Sparsity = sparsityMatrix(lftfun.LTI.D23);
D31_Sparsity = sparsityMatrix(lftfun.LTI.D31);
D32_Sparsity = sparsityMatrix(lftfun.LTI.D32);
D33_Sparsity = sparsityMatrix(lftfun.LTI.D33);
dThetadOmega_Sparsity = true(size(testo_jacobiano));

[row_zero_elements, col_zero_elements] = find(testo_jacobiano==0);
if ~isempty(find(testo_jacobiano==0))
    for ind_coord = 1:size(row_zero_elements,1)
        dThetadOmega_Sparsity(row_zero_elements(ind_coord),col_zero_elements(ind_coord))=false;
    end;
end;

VCT8_Sparsity = rowColProductBoolean(invBoolean([((diag(true(size(lftfun.LTI.C1,1),1)))|D11_Sparsity), rowColProductBoolean(D12_Sparsity,dThetadOmega_Sparsity)
                                                 D21_Sparsity,                              ((diag(true(size(lftfun.LTI.C2,1),1)))|rowColProductBoolean(D22_Sparsity,dThetadOmega_Sparsity))]),[C1_Sparsity, D13_Sparsity, D11_Sparsity
                                                                                                                                                                                     C2_Sparsity, D23_Sparsity, D21_Sparsity]);                                                 
VCT9_Sparsity = ([A_Sparsity,  B1_Sparsity
                  C3_Sparsity, D31_Sparsity]|rowColProductBoolean([B1_Sparsity,  rowColProductBoolean(B2_Sparsity,dThetadOmega_Sparsity)
                                                                   D31_Sparsity, rowColProductBoolean(D32_Sparsity,dThetadOmega_Sparsity)],[VCT8_Sparsity(:,1:size(lftfun.LTI.A,2)),VCT8_Sparsity(:,size(lftfun.LTI.A,2)+size(lftfun.LTI.D13,2)+1:end)]));

VCT10_Sparsity = ([A_Sparsity, B3_Sparsity
                  C3_Sparsity, D33_Sparsity]|rowColProductBoolean([B1_Sparsity,  rowColProductBoolean(B2_Sparsity,dThetadOmega_Sparsity)
                                                                   D31_Sparsity, rowColProductBoolean(D32_Sparsity,dThetadOmega_Sparsity)],VCT8_Sparsity(:,1:size(lftfun.LTI.A,2)+size(lftfun.LTI.D13,2))));
                                      

SparsityInformations = struct();


[row,col] = find(VCT9_Sparsity(1:size(lftfun.LTI.A,1),1:size(lftfun.LTI.A,2)));
SparsityInformations.VCT9{1,1}.NumNonZeroElements = length(row);
SparsityInformations.VCT9{1,1}.row = row;
SparsityInformations.VCT9{1,1}.col = col;

[row,col] = find(VCT9_Sparsity(size(lftfun.LTI.A,1)+1:end,1:size(lftfun.LTI.A,2)));
SparsityInformations.VCT9{2,1}.NumNonZeroElements = length(row);
SparsityInformations.VCT9{2,1}.row = row;
SparsityInformations.VCT9{2,1}.col = col;

%ind_stop_z = 0;
for n = 2 : size(lftfun.DeltaSym,1)
    i = n-1;
    ind_start_z = cell2mat(lftfun.DeltaSym(n,2));
    ind_stop_z = cell2mat(lftfun.DeltaSym(n,3));
%for i=1:size(list_TypeDelta,1)
    %ind_start_z = ind_stop_z+1;
    %ind_stop_z = ind_start_z+str2double(list_TypeDelta{i,2})-1;
    %list_TypeDelta{i,3} = ind_start_z;
    %list_TypeDelta{i,4} = ind_stop_z;

    [row,col] = find(VCT9_Sparsity(1:size(lftfun.LTI.A,1),size(lftfun.LTI.A,2)+ind_start_z:size(lftfun.LTI.A,2)+ind_stop_z));
    SparsityInformations.VCT9{1,1+i}.NumNonZeroElements = length(row);
    SparsityInformations.VCT9{1,1+i}.row = row;
    SparsityInformations.VCT9{1,1+i}.col = col;
    
    [row,col] = find(VCT9_Sparsity(size(lftfun.LTI.A,1)+1:end,size(lftfun.LTI.A,2)+ind_start_z:size(lftfun.LTI.A,2)+ind_stop_z));
    SparsityInformations.VCT9{2,1+i}.NumNonZeroElements = length(row);
    SparsityInformations.VCT9{2,1+i}.row = row;
    SparsityInformations.VCT9{2,1+i}.col = col;

end;

[row,col] = find(VCT10_Sparsity(1:size(lftfun.LTI.A,1),1:size(lftfun.LTI.A,2)));
SparsityInformations.VCT10{1,1}.NumNonZeroElements = length(row);
SparsityInformations.VCT10{1,1}.row = row;
SparsityInformations.VCT10{1,1}.col = col;

[row,col] = find(VCT10_Sparsity(size(lftfun.LTI.A,1)+1:end,1:size(lftfun.LTI.A,2)));
SparsityInformations.VCT10{2,1}.NumNonZeroElements = length(row);
SparsityInformations.VCT10{2,1}.row = row;
SparsityInformations.VCT10{2,1}.col = col;


for i=1:size(lftfun.LTI.B3,2)
    [row,col] = find(VCT10_Sparsity(1:size(lftfun.LTI.A,1),size(lftfun.LTI.A,2)+i));
    SparsityInformations.VCT10{1,1+i}.NumNonZeroElements = length(row);
    SparsityInformations.VCT10{1,1+i}.row = row;
    SparsityInformations.VCT10{1,1+i}.col = col;

    [row,col] = find(VCT10_Sparsity(size(lftfun.LTI.A,1)+1:end,size(lftfun.LTI.A,2)+i));
    SparsityInformations.VCT10{2,1+i}.NumNonZeroElements = length(row);
    SparsityInformations.VCT10{2,1+i}.row = row;
    SparsityInformations.VCT10{2,1+i}.col = col;
end

lftfun.Sparsity = SparsityInformations;


end
