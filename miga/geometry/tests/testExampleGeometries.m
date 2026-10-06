classdef testExampleGeometries <aMigaTestCase

    methods (TestParameterDefinition,Static)
        function geometryName = initializeProperty()
            geometryName = getExampleGeometries();
        end
    end

    properties(TestParameter)
        geometryName
    end

    methods(Test)
        function checkCreation(testCase,geometryName)
            fh = str2func(geometryName);
            testCase.verifyWarningFree(fh);
            geo=fh();
            testCase.verifyInstanceOf(geo,'aBasis');
        end

        function checkPlot(testCase,geometryName)
            import matlab.unittest.diagnostics.Diagnostic
            import matlab.unittest.diagnostics.FigureDiagnostic;
            fh = str2func(geometryName);
            geo=fh();
            f=testCase.invisibleFigure();
            testCase.verifyWarningFree(@()plotGeometry(geo));
            testCase.log(3,FigureDiagnostic(f));

        end       
    end
end

function out=getExampleGeometries()
[mydir,~,~] = fileparts(mfilename('fullpath'));
geodir = fullfile(mydir,'..','examples','*.m'); 
list=dir(geodir);
[~,out]=arrayfun(@(s)fileparts(s.name),list,"UniformOutput",false);
end