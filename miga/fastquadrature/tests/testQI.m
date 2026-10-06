classdef testQI < aMigaTestCase



    properties(TestParameter)
        qiBasis={bSpline(6,[zeros(1,6),0:30,30*ones(1,6)]),bSpline(2,[0,0,0,1,2,3,4,5,6]),bSpline(1,[0,0,1,2,3,4,5,6,6])}
    end



methods(Test,TestTags=["qi","bSpline","univariate"])
    function testQIisProjection(testCase,qiBasis)
        qi=qiLocal(qiBasis,qiLocalNodes(qiBasis));
        vals=qiBasis.evaluate(qi.nodes);
        testCase.verifyEqual( qi.mat*vals-eye(qiBasis.dimension),zeros(qiBasis.dimension), "AbsTol",1e-14);
    end    

    function testQIGreville(testCase,qiBasis)
        qi=qiGreville(qiBasis);
        vals=qiBasis.evaluate(qi.nodes(:)');
        testCase.verifyEqual( qi.mat*vals-eye(qiBasis.dimension),zeros(qiBasis.dimension), "AbsTol",1e-14);
    end
end

end
