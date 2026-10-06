classdef parDers < evaluationOp

    properties
        request (:,:) double
        shape   (1,:) double
    end

    methods
        function pD = parDers(request,shape)
            arguments
                request
                shape (1,:) double = []
            end
            pD.request = request;
            pD.shape = shape;

        end

        function ders = what(obj,domainDim)
            assert(size(obj.request,1)==domainDim,'MIGA:invalid partial derivative specification','The partial derivative does not match the domain dimension');
            ders = obj.request;
        end

        function pD = evaluate(obj,v,~,rangeSize)
            assert(isempty(obj.shape) || prod(obj.shape)==(numel(v)*prod(size(v{1},1:2))),...
                'parDers:InputMismatch',...
                'The provided shape does not match the data.');
            if isempty(obj.shape)
                if numel(v)<=1
                    newShape = rangeSize;
                else
                    newShape = [rangeSize,numel(v)];
                end
            else
                newShape = obj.shape;
            end
            pD = reshape(cat(2,v{:}),[newShape,size(v{1},3:4)]);
        end
    end
end