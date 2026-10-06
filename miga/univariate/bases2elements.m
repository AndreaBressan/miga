function [elements] = bases2elements(bases)
%BASES2ELEMENTS computes the coarsest subdivision of the domain in elements
% of a collection of B-spline sets and or polynomials
arguments(Repeating)
    bases (:,1) aBasis {}
end
bases=cat(1,bases{:});
assert(all(arrayfun(@(b)size(b.domain,1),bases)==1,'all'),"only works for univariate bases");
allEnds=cat(1,bases.domain);
hasKnots=arrayfun(@(b)isprop(b,'knots')|ismethod(b,'knots'),bases);
allBreak=unique(cat(1,allEnds(:),bases(hasKnots).knots));
elements=[allBreak(1:end-1),allBreak(2:end)];
end
