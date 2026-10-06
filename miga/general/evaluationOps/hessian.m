classdef hessian < evaluationOp
    methods (Static)
        function ders = what(domainDim)
            S=eye(domainDim);
            ders=S+reshape(S,domainDim,1,domainDim);
            ders=reshape(ders,domainDim,domainDim^2);
        end

        function hess = evaluate(v,domainDim,rangeSize)
            hess=reshape(cat(2,v{:}),[rangeSize,domainDim,domainDim,size(v{1},[3,4])]);
        end
    end
end