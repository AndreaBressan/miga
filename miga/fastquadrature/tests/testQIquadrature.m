classdef testQIquadrature < aIntegrationTestInterface

properties
    qiBasis
    weight
    weight2
    test
    trial
    Q
end

methods(TestClassSetup)
    function setupClass(testCase)
        testCase.qiBasis=bSpline(2,[0,0,0,1,2,3,4,5,5,5]);
        testCase.weight=subBasis(testCase.qiBasis,(1:testCase.qiBasis.dimension)');
        testCase.weight2=subBasis(testCase.qiBasis,...
            reshape(1:2*testCase.qiBasis.dimension,testCase.qiBasis.dimension,2),...
            []);
        testCase.trial=bSpline(3,[0,0,0,1,2,3,4,5]);
        testCase.test =bSpline(1,[0,0,1,3,4,5]);
        elements=bases2elements(testCase.trial,testCase.test,testCase.qiBasis);
        refQ=uElementQuadrature(elements,gaussQuad(3));
        testCase.Q = qiQuadrature(qiLocal(testCase.qiBasis,qiLocalNodes(testCase.qiBasis)),refQ);
    end
end


methods(Test,TestTags=["qi","univariate","integration"])

    function testMass(testCase)
        M=testCase.Q.integrate(testCase.weight,testCase.test,zeros(1,1,2),testCase.trial);
        Mref=testCase.Q.reference.integrate(testCase.weight,testCase.test,zeros(1,1,2),testCase.trial);
        testCase.verifyEqual(M,Mref,"AbsTol",testCase.absoluteTolerance,"RelTol",testCase.relativeTolerance);
    end

    function testWithPartialDerivative(testCase)
        pairs=reshape([0,0,1,0],[1,2,2]);
        M=testCase.Q.integrate(testCase.weight2,testCase.test,pairs,testCase.trial);
        Mref=testCase.Q.reference.integrate(testCase.weight2,testCase.test,pairs,testCase.trial);
        testCase.verifyEqual(M,Mref,"AbsTol",testCase.absoluteTolerance,"RelTol",testCase.relativeTolerance);
    end

    function testRHS(testCase)
        M=testCase.Q.integrate(testCase.weight,testCase.test);
        Mref=testCase.Q.reference.integrate(testCase.weight,testCase.test);
        testCase.verifyEqual(M,Mref,"AbsTol",testCase.absoluteTolerance,"RelTol",testCase.relativeTolerance);
    end


    function testRHSwithDerivative(testCase)
        ders=reshape([0,1],[1,2]);
        M=testCase.Q.integrate(testCase.weight2,testCase.test,ders);
        Mref=testCase.Q.reference.integrate(testCase.weight2,testCase.test,ders);
        testCase.verifyEqual(M,Mref,"AbsTol",testCase.absoluteTolerance,"RelTol",testCase.relativeTolerance);
    end

    
end

end