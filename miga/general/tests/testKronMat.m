classdef testKronMat < aMigaTestCase

methods(Test,TestTags=["utils","linearAlgebra","tensor","kronecker"])

function testKronMat1(testCase)
    A1 = randi(20,3,4);
    A2 = randi(20,7,6);
    A3 = randi(20,2,8);
    A = kronMat({A1,A2,A3});
    testCase.verifyEqual(shape(A),[3 4;7 6;2 8]);
    testCase.verifySize(A,[42,192]);
    testCase.verifyEqual(shape(A'),[4 3;6 7;8 2]);
    testCase.verifySize(A',[192,42]);
    testCase.verifyEqual(toMat(A),kron(A1,kron(A2,A3)),"AbsTol",testCase.absoluteTolerance);
end

function testKronMat2(testCase)
    B1 = triu(magic(3));
    B2 = tril(magic(4));
    B3 = triu(magic(5));
    B4 = tril(magic(6));
    B5 = triu(magic(7));
    B = kronMat(B1,B2,B3,B4,B5);
    x = [1:factorial(7)/2]';
    testCase.verifyEqual(B*x,toMat(B)*x,"AbsTol",testCase.absoluteTolerance);
    testCase.verifyEqual(B\x,toMat(B)\x,"AbsTol",2e-12);
    testCase.verifyEqual(B.transpose()*x,full(B)'*x,"AbsTol",testCase.absoluteTolerance);
    testCase.verifyEqual(B.transpose()\x,full(B)'\x,"AbsTol",testCase.absoluteTolerance);
end
end
end