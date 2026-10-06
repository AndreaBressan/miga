classdef testBSpline < aBasisTestInterface
    % testSimpleBSpline 
    %
    % Usage:
    %   runtests('testSimpleBSpline')
    %   runtests('tests/testSimpleBSpline.m')
    properties (ClassSetupParameter)
        Degree = {1, 2, 3}
        Type = {'periodic', 'open', 'open-right', 'open-left', 'open-none'}
    end
    
    methods (TestClassSetup)
        function classSetup(testCase,Degree,Type)
            switch Type
                case 'periodic'
                    periodic=true;
                    knots  = 0:4;
                    domain = [-1,4];
                case 'open'
                    periodic=false;
                    knots  = [-1*ones(1,Degree),-1:4,4*ones(1,Degree)];
                    domain = [-1,4];                    
                case 'open-right'
                    periodic=false;
                    knots  = [-1*ones(1,Degree-1),-1:4,4*ones(1,Degree)];
                    domain = [-1,4];
                case 'open-left'
                    periodic=false;
                    knots  = [-1*ones(1,Degree),-1:4,4*ones(1,max(Degree-2,0))];
                    domain = [-1,4];  
                case 'open-none'
                    periodic=false;
                    knots  = [-1*ones(1,Degree-1),-1:4,4*ones(1,Degree-1)];
                    domain = [-1,4];              end
            testCase.basis=bSpline(Degree,knots, periodic,domain);
        end
    end


    methods (Test)

        function testConstructor(testCase)
            testCase.verifyNotEmpty(testCase.basis);
        end

        function testDataReplicability(testCase)
            basis=testCase.basis;
            other=bSpline(basis.degree,basis.knots,basis.periodic,basis.domain);
            testCase.verifyEqual(basis,other);
        end

        function testPartitionOfUnity(testCase)
            basis=testCase.basis;
            xs=linspace(basis.domain(1),basis.domain(2),20);
            [v,s]=basis.evaluate(xs,{parDers(0),parDers(1),parDers(2)});
            toCheck=s.nnzRow==basis.degree+1;
            w=s.toFull(v{1});
            actual=sum(w(toCheck,:),2);
            expected=ones(sum(toCheck),1);
            testCase.verifyEqual(actual,expected,"RelTol",testCase.relativeTolerance);
            expected=zeros(sum(toCheck),1);
            w=s.toFull(v{2});
            actual=sum(w(toCheck,:),2);
            testCase.verifyEqual(actual,expected,"AbsTol",1e-15,"RelTol",testCase.relativeTolerance);
            w=s.toFull(v{3});
            actual=sum(w(toCheck,:),2);
            testCase.verifyEqual(actual,expected,"AbsTol",1e-15,"RelTol",testCase.relativeTolerance);
        end

    function testHRefinement(testCase)
        S=testCase.basis;
        newK=[0.25 0.5 0.75];
        [R, mat] = S.refine(newK);
        testCase.verifyEqual(R.dimension,S.dimension+numel(newK));
        testCase.verifyEqual(R.domain,S.domain);
        testCase.verifyEqual(R.periodic,S.periodic);
        testCase.verifySize(mat,[R.dimension,S.dimension]);

        xs=linspace(S.domain(1),S.domain(2),20);
        [vS]=S.evaluate(xs);
        [vR]=R.evaluate(xs);
        testCase.verifyEqual(vR*mat,vS,"AbsTol",testCase.absoluteTolerance,"RelTol",testCase.relativeTolerance);
    end

    function testDegreeElevation(testCase)
        S=testCase.basis;
        inc=1;
        [R, mat] = S.refine([],S.degree+inc);
        numUK=numel(unique(S.knots));
        testCase.verifyEqual(R.dimension,S.dimension+inc*(numUK+S.periodic-1));
        testCase.verifyEqual(R.domain,S.domain);
        testCase.verifyEqual(R.periodic,S.periodic);
        testCase.verifySize(mat,[R.dimension,S.dimension]);

        xs=linspace(S.domain(1),S.domain(2),20);
        [vS]=S.evaluate(xs);
        [vR]=R.evaluate(xs);
        testCase.verifyEqual(vR*mat,vS,"AbsTol",testCase.absoluteTolerance,"RelTol",testCase.relativeTolerance);
    end

    function testKRefinement(testCase)
        S=testCase.basis;
        newK=[0.25 0.5 0.75];
        numUK=numel(unique(S.knots));
        inc=1;
        [R, mat] = S.refine(newK,S.degree+inc);
        testCase.verifyEqual(R.dimension,S.dimension+inc*(numUK+S.periodic-1)+numel(newK));
        testCase.verifyEqual(R.domain,S.domain);
        testCase.verifyEqual(R.periodic,S.periodic);
        testCase.verifySize(mat,[R.dimension,S.dimension]);

        xs=linspace(S.domain(1),S.domain(2),20);
        [vS]=S.evaluate(xs);
        [vR]=R.evaluate(xs);
        testCase.verifyEqual(vR*mat,vS,"AbsTol",testCase.absoluteTolerance,"RelTol",testCase.relativeTolerance);
    end

    function testGreville(testCase)
        S=testCase.basis;
        nFun=S.dimension;
        grev=S.greville;
        knots=S.knots;
        testCase.verifySize(grev,[nFun,1]);
        width=diff(S.domain);
        if S.periodic
            knots=[knots(end-S.degree+1:end)-width ;knots];
        end
        for i=1:nFun
            gi=sum(knots(i+1:i+S.degree))/S.degree;
            if gi<S.domain(1) && S.periodic
                gi=gi+width;
            end
            testCase.verifyEqual(grev(i),gi,"AbsTol",testCase.absoluteTolerance);
        end
    end
    end
end
