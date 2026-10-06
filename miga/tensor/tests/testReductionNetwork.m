classdef testReductionNetwork < aMigaTestCase

methods(Test)
function testPlots(testCase)

dim=3;
diffusion=cat(3,repmat(eye(dim),1,dim),repelem(eye(dim),1,dim));
advection=cat(3,zeros(dim,dim),eye(dim));
reaction =zeros(dim,dim,2);
adv=cat(2,diffusion,advection,reaction);
[req,map]=tNormalizeTable(adv);
reductonplan=tReductionNetwork(req);

graph_M2O=reductonplan.graphM2O();
graph_O2M=reductonplan.graphO2M();

import matlab.unittest.diagnostics.Diagnostic
import matlab.unittest.diagnostics.FigureDiagnostic;
f=testCase.invisibleFigure();
testCase.verifyWarningFree(@()plot(graph_M2O));
testCase.log(3,FigureDiagnostic(f));
testCase.verifyWarningFree(@()plot(graph_O2M));
testCase.log(3,FigureDiagnostic(f));

end
end
end

