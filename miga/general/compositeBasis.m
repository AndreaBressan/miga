classdef compositeBasis < aBasis

    properties(GetAccess=public,SetAccess=protected)
        outerBasis    aBasis = aBasisDefault
        innerBasis    aBasis = aBasisDefault
    end

    methods

        function obj = compositeBasis(outerBasis,innerBasis)
            arguments
                outerBasis      (1,1)  aBasis
                innerBasis      (1,1)  aBasis
            end
            assert(size(outerBasis.domain,1)==prod(innerBasis.rangeSize), ...
                "The number of variables of the outer basis must be equal to the range size of the inner one.");

            obj.outerBasis = outerBasis;
            obj.innerBasis = innerBasis;
            obj.domain = innerBasis.domain;
            obj.periodic = innerBasis.periodic;
            obj.dimension = outerBasis.dimension*innerBasis.dimension;
            obj.rangeSize = outerBasis.rangeSize;
        end

        function [b]  = supports(obj)
            b = repmat(obj.innerBasis.supports(),1,1,obj.outerBasis.dimension); 
            b = permute(b,[1,2,4,3]);
            b = reshape(b,[],2,obj.dimension);
        end

    end

    methods

        function [v,a] = evaluateI(obj,pts,ders)
            [uDers,~,iOld] = unique(ders','rows','stable');
            uDers = uDers';
            ords = sum(uDers,1);
            maxOrd = max(ords);
            [~,iNew2] = sortrows([ords(:),uDers.']);
            newDers = uDers(:,iNew2);
            %maxDers = max(newDers,[],2);

            %allInnerDers = getAllDers(maxDers);
            %
            tmp = [];
            for iDer = 1:size(newDers,2)
                tmp = [tmp;getAllDers(newDers(:,iDer))'];
            end
            allInnerDers = unique(tmp,"rows")';
            %
            allOuterDers = getAllDersOrderK(size(obj.outerBasis.domain,1),maxOrd);
            [~,outerDersOrderId] = unique(sum(allOuterDers,1));
            outerDersOrderId(end+1) = size(allOuterDers,2)+1;

            if ~any(newDers(:,1))
                valueFlag = true;
                newDers(:,1) = [];
            else
                valueFlag = false;
            end
            parts = arrayfun(@(iDer)getPartitions(newDers(:,iDer)),1:size(newDers,2),'UniformOutput',false);

            [innerDers,innerSparsity] = obj.innerBasis.evaluateI(pts,allInnerDers);
            innerDers = cat(2,innerDers{:});
            innerDers = innerSparsity.toFull(innerDers);
            phyPts = reshape(innerDers(:,1,:,:),[size(innerDers,1),prod(size(innerDers,3:ndims(innerDers)))]);
            [outerDers,outerSparsity] = obj.outerBasis.evaluateI(phyPts,allOuterDers);
            outerDers = cat(2,outerDers{:});
            outerDers = reshape(outerSparsity.toFull(outerDers),[size(outerDers,1:2),innerSparsity.shape,obj.outerBasis.dimension]);
            % outerDers = permute(outerDers,[1,2,3,5,4]);
            outerDers = reshape(outerDers,[size(outerDers,1:3),obj.dimension]);
            oldGDerSize = size(innerDers,1:2);
            innerDers = reshape(innerDers,[prod(size(innerDers,1:2)),size(innerDers,3:max(ndims(innerDers),3))]);
            innerDers = repmat(innerDers,[1,1,1,obj.outerBasis.dimension]);
            innerDers = reshape(innerDers,[size(innerDers,1:2),obj.dimension]);

            v = cell(1,size(newDers,2));
            for iDer = 1:size(newDers,2)
                v{iDer} = zeros([size(outerDers,1),size(outerDers,3:ndims(outerDers))]);
                for iGamma = 1:numel(parts{iDer})
                    C = fdbCoeff(newDers(:,iDer),parts{iDer}{iGamma});
                    for iDelta = outerDersOrderId(size(parts{iDer}{iGamma},2)+1):(outerDersOrderId(size(parts{iDer}{iGamma},2)+2)-1)
                        tmp = reshape(outerDers(:,iDelta,:,:),[size(outerDers,1),size(outerDers,3:max(3,ndims(outerDers)))]);
                        comps = getComps(allOuterDers(:,iDelta));
                        for iComp = size(comps,2)
                            [~,locGDers] = ismember(parts{iDer}{iGamma}',allInnerDers','rows');
                            linId = sub2ind(oldGDerSize,comps(:,iComp),locGDers);
                            tmp = tmp .* reshape(prod(innerDers(linId,:),1),[1,size(innerDers,2:ndims(innerDers))]);
                        end
                        v{iDer} = v{iDer} + C * tmp;
                    end
                end
                v{iDer} = reshape(v{iDer},[size(v{iDer},1),1,prod(size(v{iDer},2:ndims(v{iDer})))]);
            end
            if valueFlag
                v = cat(2,{outerDers(:,1,:)},v);
            end

            iOld2 = zeros(size(iNew2));
            iOld2(iNew2) = 1:numel(iNew2);
            loc = iOld2(iOld);
            v = v(loc);
            a = fSparsity([size(v{1},3)/obj.dimension,obj.dimension]);
        end
    end

end
