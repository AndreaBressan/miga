function [sQuad] = tSideQuadrature(quad,domain,side)
%TSIDEQUADRATURE.M gives a quadrature for a side of an hypercube
arguments(Input)
    quad (1,1) tQuadrature
    domain (:,2) double
    side   (:,1) double {checkUniformDim(side,domain)}
end
arguments(Output)
    sQuad (1,1) tQuadrature
end
comps=quad.components;
pos=tPosition(domain,side);
comps(side~=0)=arrayfun(@(p)uSinglePointQuad(p),pos);
sQuad=tQuadrature(comps);
end

function b=checkUniformDim(side,domain)
b=size(side,1)==size(domain,1);
end 