classdef testCompositeBasis < aBasisTestInterface


methods(TestClassSetup)
    function setupClass(testCase)
        iB = exprBasis([0 1;0 1],["(cos(pi/4)*x1+sin(pi/4)*x2)","(-sin(pi/4)*x1+cos(pi/4)*x2)"],2);
        oB = exprBasis([0 1;0 1],["(cos(-pi/4)*x1+sin(-pi/4)*x2)","(-sin(-pi/4)*x1+cos(-pi/4)*x2)"],2);
        testCase.basis=compositeBasis(oB,iB);
    end
end

methods(Test,TestTags=["compositeBasis"])
    function testValues(testCase)
        cB=testCase.basis;
        [cBval] = cB.evaluate({[0:0.1:1],[0:0.5:1]},parDers([0 1;0 0]));
        [cubeVal] = cube(2).evaluate({[0:0.1:1],[0:0.5:1]},parDers([0 1;0 0]));
        testCase.verifyEqual(cBval(:),cubeVal(:),"AbsTol",testCase.absoluteTolerance);
    end

    function testValues2(testCase)
        iB = exprBasis([0 1;0 1],["sin(x1+x2)","sin(x1*x2)","sin(x1*x2)"],1);
        oB = exprBasis([-1 1],["exp(x1)","exp(x1)"],1);
        cBex = exprBasis([-1 1;-1 1],["exp(sin(x1+x2))","exp(sin(x1*x2))","exp(sin(x1*x2))","exp(sin(x1+x2))","exp(sin(x1*x2))","exp(sin(x1*x2))"],1);
        cB = compositeBasis(oB,iB);
        [cBval] = cB.evaluate({[0:0.1:1],[0:0.5:1]},parDers([5 2;4 3]));
        [cBexval] = cBex.evaluate({[0:0.1:1],[0:0.5:1]},parDers([5 2;4 3]));
        testCase.verifyEqual(cBval(:),cBexval(:),"AbsTol",2e-11);
    end
end
end
