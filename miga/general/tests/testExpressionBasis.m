classdef testExpressionBasis < aBasisTestInterface

methods(TestClassSetup)
    function classSetup(testCase)
        % trigonometric parametrization of a circular sector lifted in the z direction
        B=exprBasis([0,1;1,2],["sin(k1*x1)*x2","cos(k1*x1)*x2","x1+x2","k2"],4);

        testCase.basis=B.setParameters(["k1","k2"],[5,3]);
    end
end

methods(Test,TestTags=["expressionBasis"])
    function verifyValues(testCase)
        th=linspace(0,1,100);
        ra=linspace(1,2,100);
        V=testCase.basis.evaluate({th,ra});
        testCase.verifyEqual( reshape(sqrt(V(1,:).^2+V(2,:).^2),100,100)-ra,zeros(100,100),"AbsTol",testCase.absoluteTolerance);
        testCase.verifyEqual( reshape(V(4,:),100,100)-3*ones(100,100),zeros(100,100),"AbsTol",testCase.absoluteTolerance);
        V=testCase.basis.evaluate({th,ra},parDers([1,0;0,1]));
        areaWeight = reshape(prod(pageeig(V(1:2,1:2,:),'vector'),1),100,100);
        testCase.verifyEqual( areaWeight-ra*testCase.basis.parametersValue(1),zeros(100,100), "AbsTol",1e-14);
    end
end


end