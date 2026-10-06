classdef evaluationOp < matlab.mixin.Heterogeneous

    methods(Abstract)
        ders = what(domainDim)
        d = evaluate(v,domainDim,rangeSize)
    end

end