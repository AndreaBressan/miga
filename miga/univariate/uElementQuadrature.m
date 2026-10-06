classdef uElementQuadrature < aIntegration
    % UQUADRATURE
    % represents a composite quadrature on given elements
    % computed using a specified model

    properties(GetAccess=public,SetAccess=protected)
        elements (:,2) double
        models   (:,1) struct
    end

    methods
        function obj = uElementQuadrature(elements,models)
            % uElementQuadrature
            % implements a univariate composite quadrature where
            % for given integration elements and model[s].
            % In each element the given quadrature is used.
            arguments
                elements (:,2) double {}
                models   (:,1) struct {}
            end
            assert(all(elements(:,2)>=elements(:,1)),'Wrongly oriented elements in univariate quadrature');
            obj.domain=[min(elements(:)), max(elements(:))];
            obj.elements=elements;
            obj.models=models;
            if isscalar(models)
                obj.nodes=reshape( (1-models.nodes)/2.*elements(:,1)+(1+models.nodes)/2.*elements(:,2),1,[]);
                obj.weights=reshape(models.weights.*(elements(:,2)-elements(:,1))/2,1,[]);
            else
                nel=size(elements,1);
                assert(nel==numel(models),'Different number of elements and quadratures provided');
                NC=arrayfun(@(i)(-1+models(i).nodes)/2.*elements(i,1)+(1-models(i).nodes)/2.*elements(i,2),1:nel,'UniformOutput',false);
                WC=arrayfun(@(i)models(i).weights.*(elements(i,2)-elements(i,1)) ,'UniformOutput',false);
                obj.nodes=cat(2,NC{:});
                obj.weights=cat(2,WC{:});
            end
            [obj.nodes,o]=sort(obj.nodes);
            obj.weights=reshape(obj.weights(o),[],1);
        end
    end
end

