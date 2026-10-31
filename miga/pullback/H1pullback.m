classdef H1pullback < aBasis
    % H1pullback is a basis representing the pullback of an H^1 basis
    % through a geometric mapping, supporting high-order partial derivatives.
    properties(GetAccess=public,SetAccess=private)
        geo             aBasis
        coefs           aBasis
        oldRequest      (:,:,:) double
        newRequest      (:,:,:) double
        side            (:,1)   double
        dependencyMap   (1,1)   struct
    end

    methods
        function obj = H1pullback(geo,coefs,oldRequest,side)
            arguments
                geo        (1,1)   aBasis
                coefs      (1,1)   aBasis = constantFunction(1,repmat([-inf,inf],geo.rangeSize,1))
                oldRequest (:,:,:) double = double.empty
                side       (:,1)   double = zeros(size(geo.domain,1),1)
            end
            obj@aBasis();
            obj.geo = geo;
            obj.coefs = coefs;
            obj.domain = geo.domain;
            obj.periodic = geo.periodic;
            if isempty(oldRequest)
                obj.oldRequest = zeros(geo.rangeSize,1,2);
            else
                obj.oldRequest = oldRequest;
            end
            [obj.newRequest,obj.dependencyMap] = obj.computeRequest(obj.oldRequest);
            obj.dimension = size(obj.newRequest,2);
            obj.side = side;
        end

        function [b]  = supports(obj)
            b = repmat(obj.geo.domain,1,1,obj.dimension);
        end

        function [newRequest,map]=computeRequest(obj,oldRequest)
            numOldRequest=size(oldRequest,2);
            dim=size(obj.geo.domain,1); % thisize(oldRequest,1);
            order=sum(oldRequest,1);

            maxDTest  = max(order(:,:,1));
            numDTest  = nchoosek(dim+maxDTest,dim);
            maxDTrial = max(order(:,:,2));
            numDTrial = nchoosek(dim+maxDTrial,dim);

            map.max=[maxDTest,maxDTrial];
            map.num=[numDTest,numDTrial];
            map.which=false([map.num,numOldRequest]);

            for req=1:numOldRequest
                map.which(:,:,req)=order2Dep(order(1,req,1),numDTest,dim)'.*order2Dep(order(1,req,2),numDTrial,dim);
            end
            allRequired=any(map.which,3);
            map.destination(allRequired)=1:sum(allRequired,'all');

            [r,c]=find(allRequired);
            map.ders=allDers(dim,max(map.max));
            newRequest=cat(3,map.ders(:,r),map.ders(:,c));
        end

    end

    methods
        function [val,a] = evaluateI(obj,pts,ders)
            arguments
                obj    H1pullback
                pts
                ders   (:,:) double
            end
            domDim = size(obj.geo.domain,1);
            ranDim = obj.geo.rangeSize;
            nPts = npts(pts);
            assert(isequal(ders,zeros(domDim,1)),"Pullback derivatives are not implemented because derivatives of coefficients are not needed by assembly algorithms");
            
            % Generate the master list of derivatives up to order maxDerOrder
            maxDerOrder = max(obj.dependencyMap.max);
            allReqPhys = allDers(ranDim, maxDerOrder);
            Nphys = size(allReqPhys, 2);

            % Evaluate geometric map derivatives up to max required mapping derivatives
            allReqPar = obj.dependencyMap.ders;
            Npar = size(allReqPar, 2);
            if maxDerOrder>0
                geoDers = obj.geo.evaluate(pts, parDers(allReqPar,[ranDim,Npar]));
            else
                geoDers = obj.geo.evaluate(pts, parDers( allDers(domDim,1),[ranDim,domDim+1]));
            end            

            % Initialize the Transformation Matrix P (d_parametric = P * d_physical)
            P = zeros(Npar, Nphys, nPts);
            P(1, 1, :) = 1; % 0th-order mapping is identity (values map to values)
            % Construct P using Faà di Bruno's formula
            for iRow = 2:Npar
                alpha = allReqPar(:, iRow);
                parts = getPartitions(alpha);
                for iGamma = 1:numel(parts)
                    gamma = parts{iGamma};
                    k = size(gamma, 2);
                    C = fdbCoeff(alpha, gamma);

                    % Locate where the required map derivatives are stored in geoDers
                    [~, locGDers] = ismember(gamma', allReqPar', 'rows');

                    for iCol = 1:Nphys
                        beta = allReqPhys(:, iCol);
                        % Only map to physical derivatives of total order exactly equal to number of partitions (k)
                        if sum(beta) ~= k
                            continue;
                        end

                        comps = getComps(beta);
                        term = zeros(1, 1, nPts);

                        for iComp = 1:size(comps, 2)
                            comp = comps(:, iComp);
                            prod_term = ones(1, 1, nPts);
                            for j = 1:k
                                % comp(j) specifies the physical coordinate dimension (x, y, z...)
                                % locGDers(j) specifies the parametric derivative taken on that coordinate
                                prod_term = prod_term .* geoDers(comp(j), locGDers(j), :);
                            end
                            term = term + prod_term;
                        end
                        P(iRow, iCol, :) = P(iRow, iCol, :) + C * term;
                    end
                end
            end

            % Extract the 1st-order block (Jacobian transpose) to compute volume element
            Jcols=2:domDim+1; Jcols=Jcols(~obj.side);
            J=geoDers(:,Jcols,:);
            detJ=sqrt(prod(pageeig(pagemtimes(J,'transpose',J,'none')),1));

            % Vectorized inversion of P across all quadrature points
            invP = permute(pagemldivide(P,eye(ranDim+1)),[3,1,2]);

            % Evaluate basis coefficients at mapped physical points
            oldCoefsEval = reshape(obj.coefs.evaluate((geoDers(:,1,:))), nPts, []);
            val = zeros([nPts, size(obj.newRequest,2)]);

            dest=reshape(obj.dependencyMap.destination,obj.dependencyMap.num);
            for req=1:size(obj.oldRequest,2)
                indTest =findDerPos(allReqPhys,obj.oldRequest(:,req,1));
                indTrial=findDerPos(allReqPhys,obj.oldRequest(:,req,2));
                [rows,cols]=find(obj.dependencyMap.which(:,:,req));
                for e=1:numel(rows)
                    r=rows(e); c=cols(e); d=dest(r,c);
                    val(:,d)=val(:,d)+invP(:,indTest,r).*invP(:,indTrial,c).*oldCoefsEval(:,req);
                end
            end
            val = {reshape(detJ(:) .* val, 1, 1, [])};
            if nargout==2; a=fSparsity([nPts,obj.dimension]); end
        end

    end

end

function dep=order2Dep(order,num,dim)
    dep=false(1,num);
    if order>0
        last =nchoosek(order+dim,dim);
        dep(2:last)=true;
    else
        dep(1)=true;
    end
end

function loc=findDerPos(list,der)
    loc=find(all(list==der,1),1);
end

function ders = allDers(dim,maxDerOrder)
    ranges = repmat({0:maxDerOrder}, 1, dim);
    [grids{1:dim}] = ndgrid(ranges{:});
    allDers = cell2mat(cellfun(@(x) x(:).', grids, 'UniformOutput', false).');
    totalOrder = sum(allDers, 1);
    allDers = allDers(:, totalOrder <= maxDerOrder);
    [~, idx] = sort(sum(allDers, 1));
    ders = allDers(:, idx);
end