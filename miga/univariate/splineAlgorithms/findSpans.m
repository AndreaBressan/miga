function span_idx=findSpans(knots,pts,rule)
% FINDSPANS returns the knot span indeces of the xs
% -knots    is the knot vector as a column 
% -pts      is the list of points, can have only one row,
%           but as many columns or ndim-pages
% -rule     chooses the span assigned to knots 
%           rule>0 corresponds to spans of the form [a,b[, i.e. to
%                   right-continuous spline derivatives
%           rule<= corresponds to spans of the form ]a,b], i.e. to
%                   left-continuous spline derivatives
% it is possible to compute the knot span of many combinations of knots
% and xs at the same time by passing them as multidimensional arrays
% with separate sets of non-unitary sizes

arguments
    knots (:,1) double {}
    pts   (1,:) double {}
    rule  (1,1) double {} =1
end
if rule>0
    span_idx=sum(pts>=knots,1);
else
    span_idx=sum(pts> knots,1);
end
end
