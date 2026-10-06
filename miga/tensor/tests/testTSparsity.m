classdef testTSparsity < aSparsityTestInterface

methods(TestClassSetup)

    function setupClass(testCase)
        shape1 = [6,7];
        cols1  = [1,2,2,4,5,6];
        rows1  = [2,3,5,4,5,6];
        data1  = 1:6;
        M1=sparse(rows1(:),cols1(:),data1(:),shape1(1),shape1(2));
        S1=gSparsity(M1~=0);

        shape2 = [3,2];
        cols2  = [1,1,2,2];
        rows2  = [2,3,1,3];
        data2  = [7,11,13,17];
        M2=sparse(rows2(:),cols2(:),data2(:),shape2(1),shape2(2));
        S2=gSparsity(M2~=0);
        testCase.sparsity=tSparsity([S1,S2]);
    end
end

methods(Test,TestTags=["sparsity","tensor"])
    function testAgreementWithKronMatrix(testCase)
        s=testCase.sparsity;
        F=s.uSparsity;
        M=arrayfun(@(f) f.toFull((1:f.nnz)'),F,"UniformOutput",false);
        MT=1;
        nnz=1;
        for i=numel(M):-1:1
            MT =kron(MT,M{i});
            nnz=nnz(:)'.*(1:F(i).nnz)';
        end
        testCase.verifyEqual(s.toFull(nnz(:)),MT,'Error tSparsity should be the Kronekcer product');
    end
end

end


%function [t] = testTSparsity()
%t=true;
%
%shape1 = [6,7];
%cols1  = [1,2,2,4,5,6];
%rows1  = [2,3,5,4,5,6];
%data1  = 1:6;
%M1=sparse(rows1(:),cols1(:),data1(:),shape1(1),shape1(2));
%S1=gSparsity(M1~=0);
%
%shape2 = [3,2];
%cols2  = [1,1,2,2];
%rows2  = [2,3,1,3];
%data2  = [7,11,13,17];
%M2=sparse(rows2(:),cols2(:),data2(:),shape2(1),shape2(2));
%S2=gSparsity(M2~=0);
%
%M =kron(M2,M1);
%% SK=gSparsity(M~=0); would be the gSparsity of the Kronecker matrix
%S=tSparsity([S1,S2]);
%shape=size(M);
%
%data=data1(:).*data2(:)';
%rows=rows1(:)+(rows2(:)'-1)*shape1(1);
%cols=cols1(:)+(cols2(:)'-1)*shape1(2);
%
%t = testAssert(t,isequal(S.toFull(data(:),1),full(M)));
%t = testAssert(t,isequal(S.toFull(data(:)',2),reshape(full(M),[1,shape])));
%t = testAssert(t,isequal(S.toSparse(data(:)',2),M));
%t = testAssert(t,isequal(S.toSparse(repmat(data(:)',2,1),2),repelem(M,2,1)));
%
%nzs=[5;6];
%[r,c]=S.pos(nzs);
%t = testAssert(t,isequal(r,rows(nzs))&& isequal(c,cols(nzs)));
%[r,c]=S.pos();
%t = testAssert(t,isequal(r,rows(:))&& isequal(c,cols(:)));
%
%t = testAssert(t,isequal(data(:),S.toCompressed(full(M),1)));
%t = testAssert(t,isequal(data(:)',S.toCompressed(reshape(full(M),[1,shape]),2)));
%t = testAssert(t,isequal(S.nzPos(rows(nzs),cols(nzs)),nzs));
%
%t = testAssert(t,isequal(sum(M~=0,1)',S.nnzCol()));
%t = testAssert(t,isequal(sum(M~=0,2),S.nnzRow()));
%
%t = testAssert(t,isequal(sum(M~=0,1)',S.nnzCol(1:shape(2))));
%t = testAssert(t,isequal(sum(M~=0,2),S.nnzRow(1:shape(1))));
%
%end

