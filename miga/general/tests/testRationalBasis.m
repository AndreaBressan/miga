classdef testRationalBasis < aBasisTestInterface
    properties (ClassSetupParameter)
        DomainDim={1,2}
    end
    

    methods (TestClassSetup)
        function classSetup(testCase,DomainDim)
            switch DomainDim
            case 1
                testCase.basis = rationalBasis(bSpline(2,[0,0,0,1,1,1]),[1;sqrt(2)/2;1]);
            case 2
                S = bSpline(1,[0,0,1,1]);
                coefs = reshape([0,1;1,1],[2,2,1]);
                S(1:2) = subBasis(S,coefs);
                T = tBasis(S);
                Tcoefs = reshape([0,1,0,0;0,0,1,0]',[4,2,1]);
                T2 = subBasis(T,Tcoefs);
                testCase.basis = rationalBasis(T2,[2,1]');
            otherwise, error('Check testRationalBasis source code!');
            end
        end
    end


    methods (Test)
        function testCaseProperties(testCase)
            bas=testCase.basis;
            domainDim=size(bas.domain,1);
            switch domainDim
            case 1
                coeffs = reshape([1,0;1,1;0,1],3,1,2);
                f = subBasis(bas,coeffs);
                v = f.evaluate(linspace(0,1,100));
                testCase.verifyEqual(v(1,:).^2+v(2,:).^2,ones(1,100), 'AbsTol',testCase.absoluteTolerance );
            case 2
                f = subBasis(bas,[0;1]);
                hes=f.evaluate({linspace(0.5,1,10),linspace(0.5,1,10)},hessian);

                dfxxEx = @(x,y)-4*y./(x+2*y).^3;
                dfxyEx = @(x,y)2*(x-2*y)./(x+2*y).^3;
                [X,Y]=ndgrid(linspace(0.5,1,10),linspace(0.5,1,10));
                dfxxEx = dfxxEx(X(:),Y(:));
                dfxyEx = dfxyEx(X(:),Y(:));
                testCase.verifyEqual( reshape(dfxxEx,[],1),squeeze(hes(1,1,:)), 'AbsTol',testCase.absoluteTolerance );
                testCase.verifyEqual( reshape(dfxyEx,[],1),squeeze(hes(1,2,:)), 'AbsTol',testCase.absoluteTolerance );
            otherwise, error('Check testRationalBasis source code!');
            end
        end
    end
end