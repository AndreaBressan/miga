classdef uSinglePointQuad < aIntegration
    % UQUADRATURE
    % represents a composite quadrature on given elements
    % computed using a specified model

    methods
        function obj = uSinglePointQuad(pos)
            arguments
                pos (1,1) double {}
            end
            obj.nodes=pos;
            obj.weights=1;
            obj.domain=[pos,pos];
        end
    end
end