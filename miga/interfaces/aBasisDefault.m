classdef aBasisDefault<aBasis

methods
    function  [v,a] = evaluateI(~,~,~)  %#ok<STOUT>
        assert(false,'MIGA:uninitialized basis','an array of aBasis is by default initialized with aBasisDefault whose evaluation is defined to error');
    end
    function  [b]  = supports (~)       %#ok<STOUT>
        assert(false,'MIGA:uninitialized basis','an array of aBasis is by default initialized with aBasisDefault whose evaluation is defined to error');
    end
end

end
