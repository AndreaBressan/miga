function ptsIn = forceInside(pts,domain,tol) 
%FORCEINSIDE Force points inside a box within a prescribed tolerance.
%
%   PTSIN = forceInside(PTS,DOMAIN)
%   forces the points in PTS to lie inside the box specified by DOMAIN.
%   Points that are outside the domain are moved to the closest point on
%   the corresponding boundary of the box.
%
%   PTSIN = forceInside(PTS,DOMAIN,TOL)
%   uses TOL as the maximum allowed distance by which a point can be
%   modified. An error is raised if any point lies farther than TOL from
%   the domain.
%
%   Input arguments:
%       PTS    - Array of points to be checked. Its first dimension
%                specifies the spatial dimension, and must agree with the
%                number of rows of DOMAIN.
%
%       DOMAIN - Array of size [D,2] defining a D-dimensional box. Each row
%                contains the lower and upper bounds of the corresponding
%                coordinate.
%
%       TOL    - Nonnegative tolerance for the correction of points lying
%                outside DOMAIN. An error is raised if the distance between
%                an input point and its corrected point exceeds TOL.
%                Default: 10e-10.
%
%   Output arguments:
%       PTSIN - Array with the same size as PTS containing the corrected
%               points. Each coordinate is clipped to the corresponding
%               interval specified by DOMAIN.


if nargin<3; tol=10e-10; end
dim=size(pts,1);
assert(dim==size(domain,1),'The domain and the points must have the same dimension')
ptsIn=zeros(size(pts));
for d=1:dim
    ptsIn(d,:)=min(max(pts(d,:),domain(d,1)),domain(d,2));
    assert(max(abs(pts(d,:)-ptsIn(d,:)),[],"all")<=tol);
end
end

