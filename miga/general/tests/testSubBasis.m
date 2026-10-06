classdef testSubBasis < aBasisTestInterface

properties(ClassSetupParameter)
    DomainDim={1,2}
    RangeSize={zeros(1,0),2}
end

properties
    coefs
    rangeSize
    domainDim
end

methods(TestClassSetup)
    function setupClass(testCase,DomainDim,RangeSize)
        deg  = 2;
        nsub = 1;
        knots= [-ones(1,deg),linspace(-1,1,nsub+1),ones(1,deg)];
        uB(1:DomainDim)=bSpline(deg,knots);
        tB=tBasis(uB);

        newDimension=5;
        testCase.coefs = rand(tB.dimension,newDimension,prod(RangeSize));
        testCase.rangeSize=RangeSize;
        testCase.domainDim=DomainDim;
        testCase.basis = subBasis(tB,testCase.coefs,RangeSize);
    end
end

methods(Test,TestTags=["basis","subBasis"])
    function testRangeSizeAndDimension(testCase)
        testCase.verifyEqual(testCase.basis.rangeSize,testCase.rangeSize);
        testCase.verifyEqual(testCase.basis.dimension,5);
    end
end
end

%
%
%function [t] = testSubBasis()
%t=true;
%
%coefs = [0,1,0;0,0,1;1,0,0;1,1,1]';
%f = subBasis(bSpline(1,[0,0,0.5,1,1]),coefs);
%[v,a] = f.evaluate(linspace(0,1,5),parDers(0));
%r=a.toFull(v);
%t = testAssert(t,all(size(a.toSparse(v))==[5,4]) && all(r(:,end)==1));
%
%[v,a] = f.evaluate(linspace(0,1,5),parDers([0,1]));
%r=a.toFull(v);
%t = testAssert(t,all(size(r)==[2,5,4]) && all(r(1,:,end)==1)&& all(r(2,:,end)==0));
% 
%
%% vectorial subBasis
%coefs = [0,1,0;0,0,1;1,0,0;1,1,1]';
%coefs = repmat(coefs,1,1,2);
%B=bSpline(1,[0,0,0.5,1,1]);
%f2 = subBasis(B,coefs,2);
%[v,a] = f2.evaluate(linspace(0,1,5),parDers(0));
%r=a.toFull(v);
%t = testAssert(t,all(size(r)==[2,5,4]) && all(r(1,:,end)==1) && all(r(2,:,end)==1));
%
%
%coefs2 = repmat(coefs,1,1,3);
%f3 = subBasis(B,coefs2,[2,3]);
%[v,a] = f3.evaluate(linspace(0,1,5),parDers(0));
%r=a.toFull(v);
%t = testAssert(t,all(size(r)==[2,3,5,4]) && all(r(:,:,:,end)==1,'all') );
%
%
%[v,a] = f3.evaluate(linspace(0,1,5),parDers(1));
%r=a.toFull(v);
%t = testAssert(t,all(size(r)==[2,3,5,4]) && all(r(:,:,:,end)==0,'all') );
%
%
%% subBasis with multivariate functions
% 
%deg = 2;
%nsub = 1;
% 
%knots = [zeros(1,deg),linspace(0,1,nsub+1),ones(1,deg)];
%uB(1:2) = bSpline(deg,knots);
%B = tBasis(uB);
%
%coefs = zeros(2,9);
%coefs(2,:)=1;
%coefs(1,1)=1;
%coefs=coefs';
%f4 = subBasis(B,coefs,1);
%
%[v,a] = f4.evaluate([linspace(0,1,5);linspace(0,1,5)],parDers([0;0]));
%r=a.toFull(v);
%t = testAssert(t,all(size(r)==[1,5,2]) && all(r(1,:,2)==1,'all') );
%
%
%% vectorial subBasis with multivariate functions
%coefs = repmat(coefs,1,1,2);
%f5 = subBasis(B,coefs,2);
%
%[v,a] = f5.evaluate([linspace(0,1,5);linspace(0,1,5)],parDers([0;0]));
%r=a.toFull(v);
%t = testAssert(t,all(size(r)==[2,5,2]) && all(r(:,:,2)==1,'all') ); 
%
%end
