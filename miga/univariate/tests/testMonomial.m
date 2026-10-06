classdef testMonomial < aBasisTestInterface
    %testMonomial Tests for monomials basis class
    %
    % Tests the monomial basis (1, x, x^2, ..., x^n) and verifies:
    % - Interface compliance (via aBasisTestInterface)
    % - Correct evaluation and derivatives
    % - Domain and periodicity handling
    
    properties (ClassSetupParameter)
        Degree = {0, 1, 2, 3, 5}
        UsePeriodic = {false, true}
    end
    
    methods (TestClassSetup)
        function classSetup(testCase, Degree, UsePeriodic)
            if UsePeriodic
                domain = [-4, 4];
            else
                domain = [-2, 2];
            end
            testCase.basis = monomials(Degree, domain, UsePeriodic);
        end
    end
    
    methods (Test, TestTags={'monomial', 'polynomial'})
        
        
        function testDerivatives(testCase)
            % Test that derivatives are computed correctly
            % For monomials: d/dx(x^n) = n*x^(n-1)
            pts = [0.5, 1.5];
            deg = testCase.basis.degree;
            ders = arrayfun(@(o)parDers(o),0:3,"UniformOutput",false);
            [v, a] = testCase.basis.evaluate(pts, ders);
            
            % Check that we got a cell array of derivatives
            testCase.verifyClass(v, 'cell', 'Derivatives should be returned as cell array');
            testCase.verifyEqual(numel(v), numel(ders), ...
                'Number of derivative outputs should match number of derivative requests');
            
            % For each derivative order
            for der = (1:numel(ders))-1
                if der > deg
                    % Derivative beyond degree should be zero
                    testCase.verifyEqual(a.toFull(v{der+1}), zeros(npts(pts), deg+1), ...
                        sprintf('Derivative of order %d should be zero', der));
                end
            end
        end
        
        function testMonomialValues(testCase)
            % Test specific values of monomials
            % x^0 = 1, x^1 = x, x^2 = x^2, etc.
            x = 2.0;  % Test point
            [v] = testCase.basis.evaluate(x);
            expected = 0:testCase.basis.degree;
            expected = x .^ expected;            
            testCase.verifyEqual(v, expected, 'RelTol', 1e-13, ...
                'Monomial values incorrect');
        end
        
        function testFirstDerivative(testCase)           
            x = 1.5;
            [v, a] = testCase.basis.evaluate(x, {parDers(0),parDers(1)});
            v0_full = a.toFull(v{1});
            v1_full = a.toFull(v{2});
            
            % Expected values: f(x) = x^n, f'(x) = n*x^(n-1)
            n = 0:testCase.basis.degree;
            f_expected = x .^ n;
            fprime_expected = n .* x .^ (n - 1);
            
            testCase.verifyEqual(v0_full(1, :), f_expected, 'RelTol', 1e-13, ...
                'Function values incorrect');
            testCase.verifyEqual(v1_full(1, :), fprime_expected, 'RelTol', 1e-13, ...
                'First derivative values incorrect');
        end
        

        
        function testPeriodicBehavior(testCase)
            if ~testCase.basis.periodic
                testCase.verifyTrue(true, 'Skip non-periodic test');
                return;
            end

            domain = testCase.basis.domain;
            width = domain(2) - domain(1);
            x_in = domain(1) + 0.3;
            x_out = x_in + width * 2;
            
            [v_in, ~]  = testCase.basis.evaluate(x_in);
            [v_out, ~] = testCase.basis.evaluate(x_out);
            
            testCase.verifyEqual(v_in, v_out, 'RelTol', 1e-12, ...
                'Periodic wrapping failed');
        end
        
    end
    
end
