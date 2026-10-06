function [pos] = tPosition(domain,side)
%TPOSITION 
pos=domain(sub2ind(size(domain),find(side~=0),1+(side(side~=0)==1) ));
end

