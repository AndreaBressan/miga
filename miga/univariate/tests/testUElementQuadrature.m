classdef testUElementQuadrature < aIntegrationTestInterface

methods(Test)

function testGaussExpressionBasis(testCase)
f=exprBasis([0,1],["1","x1","x1.^2"],3);
els=linspace(0,1,3);
els=reshape([els(1:end-1);els(2:end)]',1,2,[]);
Q=uElementQuadrature(els,gaussQuad(3));
r=Q.integrate(f);
testCase.verifyEqual(r,[1,.5,1/3]', "AbsTol",testCase.absoluteTolerance );
end

function testGaussBSplineDeg1(testCase)
deg = 1;
knots = [0 0 1 1];
S = bSpline(deg,knots);
els=bases2elements(S);
Q2=uElementQuadrature(els,gaussQuad(3));
r=Q2.integrate(S,aBasis.empty,[],aBasis.empty);
testCase.verifyEqual(r,ones(2,1)/2,"AbsTol",testCase.absoluteTolerance);
end

function testGaussBSplineDeg2_1(testCase)
deg = 1;
knots = [0 0 0.25 0.5 0.75 1 1];
S2 = bSpline(deg,knots);
els=bases2elements(S2);
Q3=uElementQuadrature(els,gaussQuad(3));
r=Q3.integrate(S2,aBasis.empty,[],aBasis.empty);
testCase.verifyEqual(r,[1/8,1/4,1/4,1/4,1/8]',"AbsTol",testCase.absoluteTolerance);
r=Q3.integrate(S2);
testCase.verifyEqual(squeeze(r),[1/8 1/4 1/4 1/4 1/8 ]',"AbsTol",testCase.absoluteTolerance); 
end


function testGaussBSplineDeg2_2(testCase)
% Test with degree 2
deg = 2;
knots = [0 0 0 1 1 1];
S = bSpline(deg,knots);
els=bases2elements(S);
Q2=uElementQuadrature(els,gaussQuad(3));
r=Q2.integrate(S,aBasis.empty,[],aBasis.empty);
testCase.verifyEqual(r,ones(3,1)/3,"AbsTol",testCase.absoluteTolerance ); 
end

function testMassBSplineDeg1(testCase)
% Test mass matrix
deg1   = 1;
knots1 = [0 0 .25 .5 .75 1 1];
S1 = bSpline(deg1,knots1);
deg2   = 2;
knots2 = [0 0 0 1 1 1];
S2 = bSpline(deg2,knots2);
els=bases2elements(S1,S2);
Q2=uElementQuadrature(els,gaussQuad(3));
M=Q2.integrate(constantFunction(1,[min(knots1),max(knots1)]),S1,zeros(1,1,2),S2);
testCase.verifyEqual(full(sum(M,1)),[1/3, 1/3, 1/3],"AbsTol",testCase.absoluteTolerance);
testCase.verifyEqual(full(sum(M,2)),[1/8; 1/4; 1/4; 1/4; 1/8],"AbsTol",testCase.absoluteTolerance);
end

% Test rhs
function testRHS(testCase)
f=exprBasis([-1,1],'x1');
deg   = 2;
S = monomials(deg,[-1,1]);
Q=uElementQuadrature([-1,1],gaussQuad(3));
testCase.verifyEqual(Q.integrate(f,S),[0;2/3;0],"AbsTol",testCase.absoluteTolerance);
end
end
end