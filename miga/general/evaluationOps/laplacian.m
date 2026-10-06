classdef laplacian < evaluationOp
    methods (Static)
        function ders = what(domainDim)
            ders = 2*eye(domainDim);
        end

        function lapl = evaluate(v,~,rangeSize)
            for i=2:numel(v)
                v{1}=v{1}+v{i};
            end
            lapl=reshape(v{1},[rangeSize,size(v{1},[3,4])]);
        end
    end
end