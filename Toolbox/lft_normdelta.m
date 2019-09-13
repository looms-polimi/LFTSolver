function [delta] = lft_normdelta(lftfun, fullvalues)
% genera un delta normalizzato a partire da un delta full-scale

dim=size(lftfun.DeltaSym);
delta = zeros(cell2mat(lftfun.DeltaSym(end,3)),1);
jmax=dim(1,1)-1;
for j = 1 : jmax
    imin=cell2mat(lftfun.DeltaSym(1+j,2));
    imax=cell2mat(lftfun.DeltaSym(1+j,3));
    for i = imin : imax
         max = cell2mat(lftfun.DeltaSym(j+1,5));
         min = cell2mat(lftfun.DeltaSym(j+1,4));
         delta(i,1) = (fullvalues(j,1) - (max + min)/2)/((max - min)/2);
    end
end
delta = diag(delta);
end
