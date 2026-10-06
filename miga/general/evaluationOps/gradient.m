classdef gradient < evaluationOp
    methods (Static)
        function ders = what(domainDim)
            ders = eye(domainDim);
        end

        function grad = evaluate(v,domainDim,rangeSize)
            grad=reshape(cat(2,v{:}),[rangeSize,domainDim,size(v{1},[3,4])]);
        end
    end
end