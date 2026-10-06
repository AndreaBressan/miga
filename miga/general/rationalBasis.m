classdef rationalBasis < aBasis
    % rationalBasis is a rational combination of a set of
    % functions: a single function, a basis or generating set.
    % Each function is scaled by an array of user-defined weights and normalized by the
    % weighted sum: R_i(x) = w_i * φ_i(x) / (sum over j of w_j * φ_j(x)).

    properties(GetAccess=public,SetAccess=protected)
        basis           aBasis=aBasisDefault
        weights (:,:,:) double=double.empty
        ratBasis        aBasis
    end

    methods

        function obj = rationalBasis(basis,weights)
            arguments
                basis        (1,1)    aBasis
                weights      (:,:,:)  double
            end
            assert( size(weights,1)==basis.dimension, ...
                "The number of weights must be equal to the dimension of the basis.");
            assert( prod(basis.rangeSize)==1, ...
                "The underlaying basis must be scalar.");
            assert( any(weights), ...
                "The weights must not be all zero.");

            obj.basis=basis;
            obj.weights=weights;
            obj.domain=basis.domain;
            obj.periodic=basis.periodic;
            obj.dimension=basis.dimension;
            obj.rangeSize=basis.rangeSize;
            obj.ratBasis = productBasis(rescaledBasis(obj.basis,weights(:)),compositeBasis(exprBasis([-inf,inf],"1/x1"),subBasis(obj.basis,obj.weights)));
        end

        function [b]  = supports(obj)
            b = obj.basis.supports();
        end

    end

    methods

        function [v,a] = evaluateI(obj,pts,ders)
            [v,a] = obj.ratBasis.evaluateI(pts,ders);
        end

    end
end