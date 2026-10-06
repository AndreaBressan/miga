classdef aMigaTestCase < matlab.unittest.TestCase
    %MIGATESTCASE Base class for all MIGA tests
    %   Provides domain-specific assertions and utilities for testing
    %   the MIGA mathematical library.
    %
    % Usage:
    %   classdef MyTest < aMigaTestCase
    %       methods (Test)
    %           function testSomething(testCase)
    %               testCase.verifyNumericallyEqual(actual, expected);
    %           end
    %       end
    %   end
    
    properties
        relativeTolerance = 1e-13
        absoluteTolerance = 1e-15
        migaDir
        testFigureDir
    end

    methods(TestClassSetup)
        function setupClass(testCase)
            [interfaceDir,~,~] = fileparts(mfilename('fullpath'));
            testCase.migaDir=fullfile(interfaceDir,'..');
            testCase.testFigureDir=fullfile(testCase.migaDir,'Documentation','testFigures');
        end
    end    

    methods
        % utility functions to make invisible figures, saving them and displaying them
        function [f]=invisibleFigure(~,varargin)
            f=figure(varargin{:},'visible','off');
        end
    end

    %methods (TestMethodSetup)
        %function setupTest(testCase)
            %SETUPTEST Called before each test method
            % Can be overridden in subclasses with:
            %   methods (TestMethodSetup)
            %       function setup(testCase)
            %           setupTest@MigaTestCase(testCase);
            %           % ... additional setup
            %       end
            %   end
        %end
    %end
    
    methods (TestMethodTeardown)
        function teardownTest(testCase)
            %TEARDOWNTEST Called after each test method
            % Clean up any created figures, files, etc.
            close all hidden
        end
    end
    
end

