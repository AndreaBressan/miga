classdef extrudedBasis < aBasis

    properties
        basis  (1,1) aBasis = aBasisDefault
        vec    (:,1) double
    end

    methods
        function obj = extrudedBasis(basis,vec)
            arguments
                basis  (1,1) aBasis
                vec    (:,1) double
            end
            obj.basis = basis;
            obj.vec = vec;
            obj.domain = [basis.domain;0 1];
            obj.periodic = basis.periodic;
            obj.dimension = basis.dimension;
            if isempty(basis.rangeSize)
                obj.rangeSize = 2;
            else
                obj.rangeSize = basis.rangeSize + 1;
            end
        end

        function [b] = supports(obj)
            b = cat(1,obj.basis.supports,[0,1]);
        end

        function [v,a] = evaluateI(obj,pts,ders)

            dim = size(ders,1);
            oldRangeDim = prod(obj.basis.rangeSize);
            nPts = npts(pts);
            isGrid=iscell(pts);
            computeValues = ~all(ders(end,:));

            % If for any ders, we are deriving w.r.t. the new variable, we
            % do not need to evaluate the original basis.
            if computeValues
                if isGrid
                    [v,a] = evaluateI(obj.basis,pts(1:end-1),ders(1:end-1,:));
                else
                    [v,a] = evaluateI(obj.basis,pts(1:end-1,:),ders(1:end-1,:));
                end
            else
                v = cell(1,size(ders,2));
            end

            for iDer = 1:size(ders,2)
                if computeValues
                    v{iDer} = a.toFull(v{iDer});
                end
                if ders(end,iDer) > 1
                    v{iDer} = zeros(oldRangeDim+1,1,nPts*obj.dimension);
                elseif ders(end,iDer) == 1
                    v{iDer} = obj.vec.*ones(dim,1,nPts*obj.dimension);
                else
                    v{iDer} = reshape(v{iDer},oldRangeDim,1,[],obj.dimension);
                    v{iDer} = repmat(v{iDer},[1,1,nPts/size(v{iDer},3),1]);
                    v{iDer} = cat(1,v{iDer},zeros(1,1,nPts,obj.dimension));
                    if ~any(ders(1:end-1,iDer))
                        if isGrid
                            v{iDer} = v{iDer} + repelem(reshape(obj.vec.*pts{end},oldRangeDim+1,1,[],1),1,1,nPts/numel(pts{end}),obj.dimension);
                            v{iDer} = reshape(v{iDer},oldRangeDim+1,1,[]);
                        else
                            v{iDer} = v{iDer} + repmat(obj.vec.*pts(end,:),1,1,1,obj.dimension);
                        end
                    end
                    v{iDer} = reshape(v{iDer},size(v{iDer},1),1,[]);
                end
            end

            a = fSparsity([size(v{1},3)/obj.dimension,obj.dimension]);
        end

        % function [v,a] = evaluateCombinationI(obj,coefs,pts,ders)
        %
        % end

    end
end