function part = side2part(n,side)
%SIDE2PART Convert side indices to the part representation.
%
%   PART = side2part(N,SIDE)
%   converts the side indices in SIDE to the corresponding part
%   representation for an N-dimensional domain.
%
%   The side indices use the convention
%
%       2*i-1  -> lower side in direction i
%       2*i    -> upper side in direction i
%       0      -> no side
%
%   The corresponding entry in PART is -1 for a lower side, +1 for an
%   upper side, and 0 when SIDE is zero.
%
%   Input arguments:
%       N    - Dimension of the domain.
%
%       SIDE - Array containing the side indices to be converted. The
%              values 0,...,2*N are interpreted according to the convention
%              described above.
%
%   Output arguments:
%       PART - Array with the same size as SIDE. Each entry contains
%              -1, 0, or 1, indicating respectively the lower side, no
%              side, or upper side associated with the corresponding entry
%              of SIDE.

sz=[n,size(side)];
part=zeros(sz);
ids=sub2ind([n,numel(side)], max(ceil(side(:)'/2),1), 1:numel(side));
part(ids)=(1-2*mod(side(:),2)).*(side(:)~=0);
end

