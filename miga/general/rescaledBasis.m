classdef rescaledBasis < aBasis

    properties
        basis   (1,1) aBasis = aBasisDefault
        weights (:,1) double
    end

    methods
        function obj = rescaledBasis(basis,weights)
            arguments
                basis   (1,1) aBasis
                weights (:,1) double
            end
            obj.basis = basis;
            obj.weights = weights;
            obj.domain = basis.domain;
            obj.periodic = basis.periodic;
            obj.dimension = basis.dimension;
            obj.rangeSize = basis.rangeSize;
        end

        function [b] = supports(obj)
            b = obj.basis.supports;
        end

        function [v,a] = evaluateI(obj,pts,ders)
            [v,a] = evaluateI(obj.basis,pts,ders);
            [~,f] = a.pos;
            t = reshape(obj.weights(f),1,1,[]);
            v = cellfun(@(w)w.*t,v,'UniformOutput',false);
        end

        function [v,a] = evaluateCombinationI(obj,coefs,pts,ders)
            [v,a] = obj.basis.evaluateCombinationI(obj,obj.weights.*coefs,pts,ders);
        end

    end

end