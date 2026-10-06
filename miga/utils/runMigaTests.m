function [result] = runMigaTests(throwOnError,coverage_enabled,filter)
% [res]=runMigaTests(false,true,"polynomial");

arguments
    throwOnError     (1,1) logical=false
    coverage_enabled (1,1) logical=false
end
arguments(Repeating)
    filter 
end
fprintf("Collect miga tests ...\n")
[utilDir,~,~] = fileparts(mfilename('fullpath'));
migaDir = fullfile(utilDir,'..');
modules=getMigaModules();
[sourceFiles,testFiles]=getMigaSourceAndTestFiles(modules);
fprintf("Found modules:\n")
display({modules.name})

import matlab.unittest.TestSuite
import matlab.unittest.selectors.HasTag

newTests=testFiles;
suite = arrayfun(@(f) TestSuite.fromFile(fullfile(f.folder,f.name)),newTests,'UniformOutput' , false);
suite = [suite{:}];
for f=1:numel(filter)
    if ischar(filter{f}) || isstring(filter{f})
        if startsWith(filter{f}, '~')
            selector = ~HasTag(filter{f}(2:end));
        else
            selector = HasTag(filter{f});
        end
        suite = suite.selectIf(selector);
    end
end
import matlab.unittest.plugins.TestRunnerPlugin
import matlab.unittest.plugins.StopOnFailuresPlugin
import matlab.unittest.plugins.CodeCoveragePlugin
import matlab.unittest.plugins.codecoverage.CoverageReport
import matlab.unittest.plugins.TestReportPlugin;

fprintf("Running %d selected tests ...\n",numel(suite))

runner = matlab.unittest.TestRunner.withTextOutput('Verbosity', 1);

if throwOnError
    runner.addPlugin(StopOnFailuresPlugin)
end
if coverage_enabled
    testResultFolder=fullfile(migaDir,"Documentation","tests");
    ensureFolder(testResultFolder);
    reportFormat = CoverageReport(testResultFolder,"MainFile","coverage.html");
    filesPath=arrayfun(@(f)string(fullfile(f.folder,f.name)),sourceFiles);
    p1 = CodeCoveragePlugin.forFile(filesPath,"Producing",reportFormat);
    resultReport=fullfile(testResultFolder,"result.html");
    p2=TestReportPlugin.producingHTML(resultReport,'Verbosity',3);
    runner.addPlugin(p1);
    runner.addPlugin(p2);
end
result = runner.run(suite);
display(result);
if coverage_enabled
    fprintf('Use ''open(%s)'' to look at the coverage report\n',fullfile(testResultFolder,"coverage.html")); 
    fprintf('Use ''open(%s)'' to look at the result report\n',resultReport); 
end

end




function ensureFolder(path)
    if ~isfolder(path)
        mkdir(path);
    end
end


function modules=getMigaModules()
    [utilDir,~,~] = fileparts(mfilename('fullpath'));
    migaDir = fullfile(utilDir,'..');
    modules = dir(migaDir); 
    modules = modules([modules.isdir]);
    excludeList=[string("."),string(".."),string("external"),string("Documentation")];
    keep=arrayfun(@(f)~ismember( string(f.name),excludeList),modules);
    modules=modules(keep);
end

function [sourceFiles,testFiles]=getMigaSourceAndTestFiles(dirs)
    allFiles=arrayfun(@(mod) dir(fullfile(mod.folder,mod.name ,'**', '*.m')),dirs,"UniformOutput",false);
    allFiles=cat(1,allFiles{:});
    testFilesFlag = arrayfun(...
    @(f)contains(string(f.folder)+string(f.name),'test','IgnoreCase',true)...
        & ~contains(string(f.folder),'interface'),allFiles);
    testFiles=allFiles(testFilesFlag);
    sourceFiles=allFiles(~testFilesFlag);
    self= arrayfun(@(f)string(f.name),testFiles)==string("runMigaTests.m");
    testFiles(self)=[];
end

