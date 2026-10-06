function side = part2side(part) 
%PART2SIDE Convert a face representation to a side index.
%
%   SIDE = part2side(PART)
%   converts the representation PART of the boundary faces of a cuboid
%   into the side-index convention used to identify its faces.
%
%   For a D-dimensional cuboid, PART is represented by a D-by-N array.
%   Each column describes one geometric entity using entries in {-1,0,1}.
%   A nonzero entry specifies the side associated with the corresponding
%   coordinate direction.
%
%   For a face, exactly one entry of each column is nonzero. The value -1
%   identifies the lower side and +1 identifies the upper side in the
%   corresponding coordinate direction.
%
%   The side-index convention is
%
%       2*i-1  -> lower side in direction i
%       2*i    -> upper side in direction i
%
%   Input arguments:
%       PART - Array whose entries must belong to {-1,0,1}. Each column
%              represents a geometric entity using the part representation.
%              Only codimension-one entities, i.e. faces, can be converted
%              to side indices. Therefore, each column must contain at most
%              one nonzero entry.
%
%   Output arguments:
%       SIDE - Array of side indices corresponding to the columns of PART.
%              A zero value is returned for entities with no nonzero part.



if ~all(part==-1 | part==0 |part==1)
    error('part is invalid');
end
nz=part(:,:)~=0;
codim=sum(nz,1);
faces= codim==1;
if ~all(codim<=1)
    error('Only faces can be converted to side indices');
end
sz=[size(part,2:ndims(part)),1];
side=zeros(sz);
part=part(:,:);
[r,c]=find(nz);
ids=sub2ind(size(part),r,c);
side(faces)=2*r-1+(part(ids)+1)/2;
end

