classdef testConstantFunction < aBasisTestInterface

properties(ClassSetupParameter)
    RangeSize={[],[2,3]}
    Domain={[-1,1],[1,3;4,9]}
    Dimension={1,3}
end

methods(TestClassSetup)
    function classSetup(testCase,RangeSize,Domain,Dimension)
        testCase.basis=constantFunction(rand(prod(RangeSize)*Dimension,1),Domain,RangeSize,Dimension);
    end
end
end
