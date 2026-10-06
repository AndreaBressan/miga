function S = kite(dim) 
%KITE Construct a B-spline basis for a kite-shaped domain.
%
%   S = kite()
%   constructs a two-dimensional B-spline basis representing a kite-shaped
%   domain.
%
%   S = kite(DIM)
%   constructs a DIM-dimensional basis representing the corresponding
%   tensor-product kite geometry.
%
%   Input arguments:
%       DIM - Dimension of the domain. DIM must be a positive integer.
%             Default: 2.
%
%   Output arguments:
%       S - Vector-valued B-spline basis representing the parametrization
%           of the DIM-dimensional kite-shaped domain.


arguments
    dim (1,1) double {mustBeGreaterThanOrEqual(dim,1),mustBeInteger} = 2
end
 
B(1:dim)=bSpline(1,[0 0 1 1]);
tpB = tBasis(B);

n_coef=2^dim;
inds=cell(dim,1);
[inds{:}]=ind2sub(2*ones(1,dim),1:n_coef);

coef=cat(1,inds{:})-1;
coef(:,end)=3;
S=subBasis(tpB,reshape(coef',2^dim,1,dim),dim);

end

