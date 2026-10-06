classdef testFSparsity < aSparsityTestInterface

methods(TestClassSetup)

    function classSetup(testCase)
        testCase.sparsity=fSparsity([3,2]);
    end
end


end

%function [t] = testFSparsity()
%t=true;
%
%S=fSparsity([3,2]);
%data=1:6;
%
%t = testAssert(t,isequal(S.toFull(data(:),1),reshape(data,3,2)));
%t = testAssert(t,isequal(S.toFull(data,2),reshape(data,1,3,2)));
%t = testAssert(t,isequal(S.toSparse(data,2),sparse(reshape(data,3,2))));
%t = testAssert(t,isequal(S.toSparse(data,2),sparse(reshape(data,3,2))));
%[r,c]=S.pos([5,6]);
%t = testAssert(t,isequal(r,[2;3])&& isequal(c,[2;2]));
%[r,c]=S.pos();
%t = testAssert(t,isequal(r,reshape(repmat((1:3)',1,2),[],1))&& isequal(c,reshape(repmat(1:2,3,1),[],1)));
%
%t = testAssert(t,isequal(data,S.toCompressed(S.toFull(data,2),2)));
%t = testAssert(t,isequal(data(:),S.toCompressed(S.toFull(data(:),1),1)));
%
%t = testAssert(t,isequal(S.nzPos([3,3],[2,1]),data([6;3])'));
%
%end

