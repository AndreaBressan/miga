classdef testADR < aMigaTestCase

methods(Test)

function testExample(testCase)
    import matlab.unittest.diagnostics.Diagnostic;
    import matlab.unittest.diagnostics.FigureDiagnostic;
    f=testCase.invisibleFigure();
    testCase.verifyWarningFree(@exampleADReq);
    testCase.log(3,FigureDiagnostic(f));
end

end
end