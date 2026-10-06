function S = cube(dim,lengthSide) 
%CUBE Construct a multilinear basis for a hypercube or parallelepiped.
%
%   S = cube()
%   constructs a three-dimensional unit cube.
%
%   S = cube(DIM)
%   constructs a DIM-dimensional unit hypercube, with side length one in
%   every coordinate direction.
%
%   S = cube(DIM,LENGTHSIDE)
%   constructs a DIM-dimensional parallelepiped whose side lengths are
%   specified by LENGTHSIDE.
%
%   Input arguments:
%       DIM        - Dimension of the cube or parallelepiped. Default: 3.
%
%       LENGTHSIDE - Side lengths of the parallelepiped. It can be either
%                    a scalar, in which case the same side length is used
%                    in every direction, or a vector with DIM entries
%                    specifying the side length in each direction.
%                    Default: ones(DIM,1).
%
%   Output arguments:
%       S - Basis object representing the multilinear parametrization of
%           the DIM-dimensional cube or parallelepiped. The basis is
%           constructed as a tensor-product basis of linear B-splines.

switch nargin
    case 0
        dim=3;
        lengthSide=ones(dim,1);
    case 1
        lengthSide=ones(dim,1);
    case 2
        if numel(lengthSide)~=dim
            if numel(lengthSide)==1
                lengthSide=ones(dim,1)*lengthSide;
            else
                error('Wrong size parameter');
            end
        end
    otherwise
            error('Wrong number of arguments');
end
            
B(1:dim)=bSpline(1,[0 0 1 1]);
tpB = tBasis(B);

nel =prod(repmat(2,1,dim)); 

C=zeros([ dim, repmat(2,1,dim)]);
idx=cell(1,dim+1);
idx(:)={':'};
for i=1:dim
    idx{1}=i;
    idx{i+1}=2;
    C(idx{:})=lengthSide(i);
    idx{i+1}=':';
end
C =reshape(C,[dim, nel]);

coef=zeros([nel,1,dim]);
for i=1:dim
    coef(1:nel,1,i)=C(i,:)';
end
 
S=subBasis(tpB,coef,dim);
end

