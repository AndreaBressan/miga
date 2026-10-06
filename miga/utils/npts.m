function [num] = npts(pts)
%NPTS returns the number of points described by the pts parameter
if iscell(pts) % cartesian grid
    num=prod(cellfun(@npts,pts));
else
    num=size(pts(:,:),2);
end
end

