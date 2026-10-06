classdef powerp < evaluationOp
    properties
        derivative
        p
    end
    methods (Static)
        function obj=powerp(p,deriv)
            arguments
                p      (1,1) double {} =2
                deriv  (:,1) double {} =[]
            end
            obj.p=p;
            obj.derivative=deriv;
        end 

        function ders = what(domainDim)
            if ~isempty(obj.derivative); ders=obj.derivative; else, ders=zeros(domainDim,1); end
        end

        function powp = evaluate(v,~,rangeSize)
            powp=reshape(abs(v{1}).^obj.p,[rangeSize,size(v{1},[3,4])]);
        end
    end
end