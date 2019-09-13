function [Mb] = sparsityMatrix(M)
%SPARSITY_MATRIX Creates the sparsity matrix of M
%   The sparsity matrix is a boolean matrix with 'TRUE' where M has element
%   not equal to zero; and 'FALSE' instead. Very useful if M is sparse!!
inds=find(M);
Mb=false(size(M));
Mb(inds)=true; %#ok<FNDSB>
end

