function [boundaries] = tboundary(boxPart)
%TBOUNDARY given the specification of a box part returns the specifications
% of its boundaries as a matrix

arguments
    boxPart (:,1) double
end

full=boxPart==0;
num=sum(full);
pos=find(full);
boundaries=repmat(boxPart,1,2*num);
for i=1:num
    boundaries(pos(i),2*(i-1)+[1,2])=[-1,1];
end
end
