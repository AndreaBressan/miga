function S = curvedKite() 
%CURVEDKITE Construct a B-spline basis for a curved kite-shaped domain.
%
%   S = curvedKite()
%   constructs a two-dimensional vector-valued B-spline basis representing
%   a curved kite-shaped domain.
%
%   Output arguments:
%       S - Vector-valued B-spline basis defining the parametrization of
%           the curved kite-shaped domain.


B(1:2) = bSpline(2,[0 0 0 1 1 1]);
tpB = tBasis(B);
C = cat(3,[0 1/2 1; 0 0 0; 0 1/2 1],[0 3/2 3; 1/2 3/2 5/2; 1/2 1 3/2],[0 5/2 5; 1 3 5; 1 3/2 2]);

coefs=C(1:3,:,:,:);
coefs=permute(coefs,[2,3,1]);
coefs=reshape(coefs,[],1,3);
S = subBasis(tpB,coefs);
end

