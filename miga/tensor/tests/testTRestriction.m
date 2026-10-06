classdef testTRestriction < aMigaTestCase

methods(Test)

function test2DExpression(testCase)
domain = [-7 1;-2 3]; 
source = exprBasis(domain,"1+x1*x2");
side = [-1; 0];  
restrictedSource = tRestrict(source,side);
testCase.verifyEqual(restrictedSource.expression,"1+-7*x1");
testCase.verifyEqual(restrictedSource.domain,[-2 3]);
end

function test3DExpression(testCase)
domain = [-7 1;-2 3; 1 2]; 
source = exprBasis(domain,"1+3*x2+x3^2-x1");
side = [-1; 0; 1];  
restrictedSource = tRestrict(source,side);
testCase.verifyEqual(restrictedSource.expression,"1+3*x1+2^2--7");
testCase.verifyEqual(restrictedSource.domain,[-2 3]);

side = [0; 0; 1];  
restrictedSource = tRestrict(source,side);
testCase.verifyEqual(restrictedSource.expression,"1+3*x2+2^2-x1");
testCase.verifyEqual(restrictedSource.domain,[-7 1;-2 3]);
end

function test2DBSpline(testCase)
dim = 2; 
deg=3;
nsub=4;
knots = [zeros(1,deg),linspace(0,1,nsub+1),ones(1,deg)];
uB(1:dim) = bSpline(deg,knots);
B = tBasis(uB);
side = [1; 0]; 
restrictedB = tRestrict(B,side);
testCase.verifyEqual(restrictedB.domain,[0 1]);
end
end
end