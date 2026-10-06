classdef productBasis < aBasis

    properties(GetAccess=public,SetAccess=protected)
        terms    (1,:) aBasis = aBasisDefault
    end

    methods

        function obj = productBasis(terms)
            arguments (Repeating)
                terms (:,1) aBasis {}
            end
            obj.terms=flip(cat(1,terms{:}));

            for iB = 2:numel(obj.terms)
                assert(isequal(obj.terms(1).domain,obj.terms(iB).domain), ...
                    "The domain of all the bases must be equal.");
                %assert(isequal(obj.terms(1).rangeSize,obj.terms(iB).rangeSize), ...
                %  "The range size of the two basis must be equal.");
            end

            obj.domain = obj.terms(1).domain;
            obj.periodic = prod([obj.terms.periodic]);
            obj.dimension = prod([obj.terms.dimension]);
            rangeDims = arrayfun(@(x)numel(x.rangeSize),obj.terms);
            rangeDim = max(rangeDims);
            rangeSizes = arrayfun(@(x) [x.rangeSize,ones(1,rangeDim-numel(x.rangeSize))],obj.terms,'UniformOutput',false);
            rangeSizes = cat(1,rangeSizes{:});
            obj.rangeSize = max(rangeSizes,[],1);
        end

        function [b]  = supports(obj)
            b = obj.terms(1).supports();
            for iT = 2:numel(obj.terms)
                tmp = obj.terms(iT).supports();
                tmp = reshape(tmp,[size(tmp,1:2),1,size(tmp,3:ndims(tmp))]);
                tmp1 = reshape(max(b(:,1,:,:),tmp(:,1,:,:)),size(b,1),1,[]);
                tmp2 = reshape(min(b(:,2,:,:),tmp(:,2,:,:)),size(b,1),1,[]);
                b = cat(2,tmp1,tmp2);
            end
        end

    end

    methods

        function [v,a] = evaluateI(obj,pts,ders)
            %% Removing duplicated derivatives and sorting
            ords = sum(ders,1);
            [~,iNew2] = sortrows([ords(:),ders.']);
            newDers = ders(:,iNew2);
            %maxDers = max(newDers,[],2);
            %% Determining needed derivatives
            tmp = [];
            for iDer = 1:size(newDers,2)
                tmp = [tmp;getAllDers(newDers(:,iDer))'];
            end
            allDers = unique(tmp,"rows")';

            %% Precomputing
            parts = cell(1,size(newDers,2));
            derId = cell(1,size(newDers,2));
            for iDer = 1:size(newDers,2)
                parts{iDer} = getSplits(newDers(:,iDer),numel(obj.terms));
                tmp = reshape(parts{iDer},size(parts{iDer},1),[]);
                [~,tmp]=ismember(tmp',allDers','rows');
                derId{iDer} = reshape(tmp,size(parts{iDer},2:3));
            end

            %% Evaluating needed derivatives
            D = numel(obj.terms);
            termsDers = cell(1,D);
            termsSparsity = cell(1,D);
            termDims = [obj.terms.dimension];

            for iT=1:D
                [termsDers{iT},termsSparsity{iT}] = obj.terms(iT).evaluateI(pts,allDers);
                termsDers{iT} = cat(2,termsDers{iT}{:});
                termsDers{iT} = termsSparsity{iT}.toFull(termsDers{iT});

                s1 = size(termsDers{iT}, 1);
                s2 = size(termsDers{iT}, 2);
                s3 = size(termsDers{iT}, 3);

                newShape = [s1, s2, s3, ones(1, iT-1), termDims(iT), ones(1, D-iT)];
                termsDers{iT} = reshape(termsDers{iT}, newShape);
            end

            %% Assembling
            v = cell(1,size(newDers,2));
            for iDer = 1:size(newDers,2)
                coeff = arrayfun(@(iS) multinomialMulti(newDers(:,iDer), parts{iDer}(:,:,iS)), 1:size(parts{iDer},3));
                coeff = reshape(coeff, 1, [], 1);
                prod_comps = termsDers{1}(:, derId{iDer}(1,:), :, :, :, :);
                for iT = 2:D
                    comp_i = termsDers{iT}(:, derId{iDer}(iT,:), :, :, :, :);
                    prod_comps = prod_comps .* comp_i;
                end
                tmp = sum(coeff .* prod_comps, 2);
                v{iDer} = reshape(tmp, size(tmp,1), 1, []);
            end
            %% Restoring original order of requested derivatives
            iOld2 = zeros(size(iNew2));
            iOld2(iNew2) = 1:numel(iNew2);
            v = v(iOld2);
            a = fSparsity([size(v{1},3)/obj.dimension,obj.dimension]); % to improve!
        end
    end

end
