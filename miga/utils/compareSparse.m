function [absErr,relErr]=compareSparse(A,B)
% compute the drowstance of the nonzero matrix entries of A or B. 
% If for each i,j: A(i,j)=0 or B(i,j)~=0 the 
% the absolute error rows given by  
%
% abs( A(i,j)-B(i,j)),
%
% and the relative error by 
%
% abs(A(i,j)-B(i,j)) / max(abs(A(i,j),B(i,j)),
%
% the maximums over i,j are then computed and assigned to the output.

arguments
    A (:,:) double {issparse}
    B (:,:) double {issparse}
end

dMat=abs(A-B);
if nnz(dMat)==0
    absErr=0;
    relErr=0;
    return;
end
if nargout==1
    absErr=max(nonzeros(dMat));
    return;
end
maxMat=max(abs(A),abs(B));
[rows,cols,aErr]=find(dMat);
begs=[1;find(cols(1:end-1)~=cols(2:end))+1;numel(rows)+1];
rDen=zeros(size(rows));
for ci=1:numel(begs)-1
    rDen(begs(ci):begs(ci+1)-1)=maxMat(rows(begs(ci):begs(ci+1)-1), cols(begs(ci)));
end
rErr=aErr./rDen;
absErr=max(aErr);
relErr=max(rErr);
end