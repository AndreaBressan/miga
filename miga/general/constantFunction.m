classdef constantFunction < aBasis
    % constantFunction represents a set of functions whose value is the
    % same in all points

    properties(GetAccess=public,SetAccess=protected)
        value
    end

    methods
        function obj = constantFunction(value,domain,rangeSize,dimension)
            arguments
                value               double
                domain     (:,2)    double
                rangeSize  (:,1)    double = size(value)
                dimension  (1,1)    double {mustBeInteger,mustBePositive} ...
                                        = numel(value)/prod(rangeSize)
            end
            assert(~isempty(value),"the argument value cannot be empty")
            assert(prod(rangeSize)*dimension==numel(value),"incompatible rangeSize and dimension are given")
            obj.rangeSize=rangeSize;
            obj.dimension=dimension;
            obj.domain=domain;
            obj.value=reshape(value,[prod(obj.rangeSize),1,1,obj.dimension]);
            obj.periodic=false(1,size(obj.domain,1));
        end

        function [b]  = supports (obj)
            b=repmat(obj.domain,[1,1,obj.dimension]);
        end
    end

    methods
        function [v,a] = evaluateI(obj,pts,what)
            nPts=npts(pts);
            isZero=max(what,1)==0;
            a=fSparsity([nPts,obj.dimension]);
            oshape=[prod(obj.rangeSize),1,nPts*obj.dimension];
            function val=repeat(z)
                if z
                    val=zeros(oshape);
                else
                    val=reshape(repmat(obj.value,[1,1,nPts,1]),oshape);
                end
            end
            v=arrayfun(@repeat,isZero,"UniformOutput",false ); 
        end
    end

end

