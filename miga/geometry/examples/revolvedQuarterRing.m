function S = revolvedQuarterRing()
%REVOLVEDQUARTERRING Construct a rational basis for a revolved quarter ring.
%
%   S = revolvedQuarterRing()
%   constructs a rational tensor-product basis representing a quarter
%   ring obtained by revolving a suitable planar spline geometry.
%
%   The basis is constructed from quadratic and linear B-spline
%   components. Rational weights are used to represent the curved geometry
%   exactly.
%
%   Output arguments:
%       S - Basis representing the parametrization of the revolved quarter
%           ring. The resulting basis is vector-valued and includes the
%           rational factor required by the curved geometry.


B(1:2)=bSpline(2,[0 0 0 1 1 1]);
B(3)=bSpline(1,[0,0,1,1]);
B(2:3)=B([3,2]);

w = [1 sqrt(2)/2 1 sqrt(2)/2 1/2 sqrt(2)/2 1 sqrt(2)/2 1 1 sqrt(2)/2 1 sqrt(2)/2 1/2 sqrt(2)/2 1 sqrt(2)/2 1];
coefs = [1,sqrt(2)/2,1,sqrt(2)/2,1/2,sqrt(2)/2,0,0,0,2,sqrt(2),2,sqrt(2),1,sqrt(2),0,0,0;...
    0,0,-1,sqrt(2)/2,1/2,-sqrt(2)/2,1,sqrt(2)/2,-1,0,0,-1,sqrt(2),1,-sqrt(2)/2,2,sqrt(2),-1;...
    0,sqrt(2)/2,1,0,1,sqrt(2),0,sqrt(2),2,0,sqrt(2)/2,1,0,3/2,3/2*sqrt(2),0,3/2*sqrt(2),3]'./w';
coefs = reshape(coefs,3,3,2,3);
coefs = permute(coefs,[1,3,2,4]);
coefs = reshape(coefs,[],1,3);

tpB = tBasis(B);
S=productBasis(subBasis(tpB,coefs.*w',3),compositeBasis(exprBasis([-inf,inf],"1/x1"),subBasis(tpB,w')));
end

