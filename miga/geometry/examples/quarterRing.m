function S = quarterRing(r1,r2,height) 
%QUARTERRING Construct a quarter-ring basis in two or three dimensions.
%
%   S = quarterRing(R1,R2)
%   constructs a basis representing a two-dimensional quarter ring with
%   inner radius R1 and outer radius R2.
%
%   S = quarterRing(R1,R2,HEIGHT)
%   constructs a three-dimensional quarter ring obtained by extending the
%   two-dimensional geometry in the third direction over a height HEIGHT.
%
%   Input arguments:
%       R1     - Inner radius of the quarter ring. Must be nonnegative.
%                Default: 1.
%
%       R2     - Outer radius of the quarter ring. Must be positive.
%                Default: 2.
%
%       HEIGHT - Height of the three-dimensional quarter ring. If omitted
%                or empty, a two-dimensional quarter ring is constructed.
%                If specified, HEIGHT must be a positive scalar.
%
%   Output arguments:
%       S - Basis representing the parametrization of the quarter ring.
%           In the two-dimensional case, S represents the planar quarter
%           annular region. In the three-dimensional case, S represents
%           its extrusion in the third coordinate direction over the
%           interval [0,HEIGHT].

arguments
    r1 (1,1) double {mustBeNonnegative}  = 1
    r2 (1,1) double {mustBePositive}     = 2
    height (:,1) double {mustBePositive} = []
end
assert(numel(height)<=1, 'MIGA:invalid parameter', "Height can only be a scalar or empty");
if isempty(height); dim=2; else, dim=3; end
   
B(1)=bSpline(1,[0,0,1,1]);
B(2)=bSpline(2,[0 0 0 1 1 1]);

if dim==3
    B(3)=bSpline(1,[0,0,1,1]);
    C =[r1,r2,r1,r2,0, 0, r1,r2,r1,r2,0, 0;
        0, 0, r1,r2,r1,r2,0, 0, r1,r2,r1,r2 ;
        0, 0, 0, 0, 0, 0, height*ones(1,6)];
    w = [1,1,sqrt(2)/2,sqrt(2)/2,1,1,...
         1,1,sqrt(2)/2,sqrt(2)/2,1,1];
else
    w = [1,1,sqrt(2)/2,sqrt(2)/2,1,1];
    C = [r1,r2,r1,r2,0, 0;
         0, 0, r1,r2,r1,r2  ];
end

tpB = tBasis(B);
% nurbs = rationalBasis(tpB,w);

C =reshape(C,dim,[]);
coefs=zeros([tpB.dimension,1,dim]);
for i=1:dim
    coefs(1:tpB.dimension,1,i)=C(i,:)';
end
%S=subBasis(nurbs,coefs,dim);
S=productBasis(subBasis(tpB,coefs.*w',dim),compositeBasis(exprBasis([-inf,inf],"1/x1"),subBasis(tpB,w')));
end

