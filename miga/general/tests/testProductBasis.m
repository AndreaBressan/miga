classdef testProductBasis < aBasisTestInterface

methods(TestClassSetup)
    function setupClass(testCase)
        b(1) =  exprBasis([0 1],["x1"]);
        b(2) = exprBasis([0 1],["x1^5"]);
        testCase.basis = productBasis(b);
    end
end

methods(Test,TestTags=["basis","productBasis"])
    function testValues1(testCase)
        pB = testCase.basis;
        pBVal = pB.evaluate(1/2,parDers([5,3,3]));
        testCase.verifyEqual(pBVal,[360;15;15],"AbsTol",testCase.absoluteTolerance);
    end
    function testValues2(testCase)
        b(1) = exprBasis([0 1],["x1","2*x1","3*x1"]);
        b(2) = exprBasis([0 1],["x1^2","2*x1^2"]);
        c = exprBasis([0 1],["x1^3","5*x1^3"]);
        pBex = exprBasis([0 1],["x1^6","5*x1^6","2*x1^6","10*x1^6","2*x1^6","10*x1^6","4*x1^6","20*x1^6","3*x1^6","15*x1^6","6*x1^6","30*x1^6"]);
        pB = productBasis(b,c);
        pBVal = pB.evaluate(0:0.1:1,parDers([5,3,3]));
        pBValEx = pBex.evaluate(0:0.1:1,parDers([5,3,3]));
        testCase.verifyEqual(pBVal(:),pBValEx(:),"AbsTol",1e-11);
    end
end
end
