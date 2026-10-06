classdef preEvaluatedFunction < aBasis
    % preEvaluatedFunction is a trick to pass computed values at the
    % the quadrature points to integrate

    properties(GetAccess=public,SetAccess=protected)
        values
    end

    methods
        function obj = preEvaluatedFunction(values,domain,rangeSize,dimension)
            arguments
                values              double
                domain     (:,2)    double
                rangeSize  (:,1)    double {mustBePositive,mustBeInteger} = []
                dimension  (1,1)    double {mustBePositive,mustBeInteger} = 1
            end
            assert(~isempty(values),"the argument value cannot be empty")
            obj.rangeSize=rangeSize;
            obj.dimension=dimension;
            obj.domain=domain;
            obj.values=reshape(values,prod(obj.rangeSize),1,[]);
            obj.periodic=false(1,size(obj.domain,1));
        end

        function [b]  = supports (obj)
            b=repmat(domain,[1,1,obj.dimension]);
        end

        function [v,a] = evaluateI(obj,pts,what)
            nPts=npts(pts);
            assert(all(what==0),"Cannot compute derivatives of a function known at quadrature points");
            assert(nPts==size(obj.values,3)/obj.dimension, "Precomputed values used on a different set of points");
            v={obj.values};
            if nargout==2; a=fSparsity([nPts,obj.dimension]);end
        end
    end

end

