classdef testRSparsity < aSparsityTestInterface

methods(TestClassSetup)

    function classSetup(testCase)
        shape = [6,7];
        cols  = [1,2,2,4,5,6];
        rows  = [2,3,5,4,5,6];
        data  = 1:6;
        M=sparse(rows(:),cols(:),data(:),shape(1),shape(2));
        testCase.sparsity=rSparsity(gSparsity(M~=0),2*(1:shape(1)),shape.*[2,1]);
    end
end


end


%function [t] = testGSparsity()
%t=true;
%
%shape = [6,7];
%cols  = [1,2,2,4,5,6];
%rows  = [2,3,5,4,5,6];
%data  = 1:6;
%M=sparse(rows(:),cols(:),data(:),shape(1),shape(2));
%S=gSparsity(M~=0);
%
%t = testAssert(t,isequal(S.toFull(data(:),1),full(M)));
%t = testAssert(t,isequal(S.toFull(data,2),reshape(full(M),[1,shape])));
%t = testAssert(t,isequal(S.toSparse(data,2),M));
%t = testAssert(t,isequal(S.toSparse(repmat(data,2,1),2),repelem(M,2,1)));
%
%nzs=[5;6];
%[r,c]=S.pos(nzs);
%t = testAssert(t,isequal(r,rows(nzs)')&& isequal(c,cols(nzs)'));
%[r,c]=S.pos();
%t = testAssert(t,isequal(r,rows(:))&& isequal(c,cols(:)));
%
%t = testAssert(t,isequal(data(:),S.toCompressed(full(M),1)));
%t = testAssert(t,isequal(data ,S.toCompressed(reshape(full(M),[1,shape]),2)));
%t = testAssert(t,isequal(S.nzPos(rows(nzs),cols(nzs)),nzs));
%
%t = testAssert(t,isequal(sum(M~=0,1)',S.nnzCol()));
%t = testAssert(t,isequal(sum(M~=0,2),S.nnzRow()));
%
%t = testAssert(t,isequal(sum(M~=0,1)',S.nnzCol(1:shape(2))));
%t = testAssert(t,isequal(sum(M~=0,2),S.nnzRow(1:shape(1))));
%
%end

