classdef testTBasis < aBasisTestInterface
    % testSimpleBSpline 
    %
    % Usage:
    %   runtests('testSimpleBSpline')
    %   runtests('tests/testSimpleBSpline.m')
    properties (ClassSetupParameter)
        Case={1,2,3,4,5,6,7,8,9,10,11,12}
    end
    
    properties
        Degree
        Periodic
        bases
        grid
    end

    methods (TestClassSetup)
        function classSetup(testCase,Case)
            if Case<=4
                testCase.Degree   = [1,2]';
                testCase.Periodic = (bitget(uint32(Case-1),[1,2])==1)';
            else
                testCase.Degree= [3,4,5]';
                testCase.Periodic = (bitget(uint32(Case-5),[1,2,3])==1)';
            end

            testCase.assumeEqual(numel(testCase.Degree), numel(testCase.Periodic), 'Skip all tests if dimension do not agree');
            testCase.bases=arrayfun(@(d,p)monomials(d,[-d,d+p],p), testCase.Degree,testCase.Periodic);
            testCase.basis=tBasis(testCase.bases);
            domain=testCase.basis.domain;
            testCase.grid=arrayfun(@(i)linspace(domain(i,1),domain(i,2),i+3 ), 1:size(domain,1),"UniformOutput",false);
        end
    end


    methods (Test)
        function testSetProperties(testCase)
            dom=cat(1,testCase.bases.domain);
            testCase.verifyEqual(testCase.basis.domain, dom);
            testCase.verifyEqual(testCase.Degree,[testCase.basis.components.degree]');
            testCase.verifyEqual(testCase.Periodic,testCase.basis.periodic);
        end
        function testDimension(testCase)
            testCase.verifyEqual(prod([testCase.basis.components.dimension]),testCase.basis.dimension);
        end
        function testLaplacianIsHessianTrace(testCase)
            S=testCase.basis;
            dim=size(S.domain,1);
            [v,~] = S.evaluate(testCase.grid,{parDers(eye(dim)),parDers(2*eye(dim)),laplacian,hessian});
            testCase.verifyEqual(squeeze(sum(v{2},1)),v{3}',"AbsTol",testCase.absoluteTolerance,"RelTol",testCase.relativeTolerance);
            tmp=0;
            for d=1:dim
                tmp=tmp+squeeze(v{4}(d,d,:));
            end
            testCase.verifyEqual(tmp,v{3},"AbsTol",testCase.absoluteTolerance,"RelTol",testCase.relativeTolerance);
        end
        function testLaplacianIsHessianTraceForFunction(testCase)
            S=testCase.basis;
            coefs=eye(S.dimension); 
            dim=size(S.domain,1);
            [v,~] = S.evaluateCombination(coefs,testCase.grid,{parDers(eye(dim)),parDers(2*eye(dim)),laplacian,hessian});
            testCase.verifyEqual(squeeze(sum(v{2},1)),v{3}',"AbsTol",testCase.absoluteTolerance,"RelTol",testCase.relativeTolerance);
            tmp=0;
            for d=1:dim
                tmp=tmp+squeeze(v{4}(d,d,:));
            end
            testCase.verifyEqual(tmp,v{3},"AbsTol",testCase.absoluteTolerance,"RelTol",testCase.relativeTolerance);
        end

        function testEvaluationOfMultipleThings(testCase)
            S=testCase.basis;
            dim=size(S.domain,1);
            [~,vp2,lap,hes,~] = S.evaluate(testCase.grid,parDers(eye(dim)),parDers(2*eye(dim)),laplacian,hessian);
            testCase.verifyEqual(squeeze(sum(vp2,1)),lap');
            tmp=0;
            for d=1:dim
                tmp=tmp+squeeze(hes(d,d,:));
            end
            testCase.verifyEqual(tmp,lap);
        end

    end
end


