classdef testLagrange < aBasisTestInterface
    %testLagrange Tests for Lagrange basis class
    %
    % Tests the Lagrange polynomial basis (dual to evaluation at nodes)
    % and verifies:
    % - Interface compliance (via aBasisTestInterface)
    % - Kronecker delta property: L_i(x_j) = delta_ij
    % - Correct derivatives
    % - Domain handling
    
    properties (ClassSetupParameter)
        NumNodes = {2, 3, 5, 7}
        UsePeriodicFlag = {false}  % Lagrange typically doesn't use periodic
    end
    
    methods (TestClassSetup)
        function classSetup(testCase, NumNodes, UsePeriodicFlag)
            % Create Lagrange basis with specified number of nodes
            % Use nodes in [-1, 1] for nice test values
            knots = linspace(-1, 1, NumNodes)';
            testCase.basis = lagrange(knots, [-1, 1], UsePeriodicFlag);
        end
    end
    
    methods (Test, TestTags={'lagrange', 'polynomial'})
        
        
        function testKroneckerDeltaProperty(testCase)
            % Most important property: L_i(x_j) = delta_ij
            % Lagrange functions equal 1 at their node, 0 at other nodes
            
            knots = testCase.basis.knots';
            [v, a] = testCase.basis.evaluate(knots);
            
            % Convert from sparse format to full
            v_full = a.toFull(v);
            
            % Check against identity matrix
            expected = eye(numel(knots));
            testCase.verifyEqual(v_full, expected, 'AbsTol', 1e-13, ...
                'Kronecker delta property violated: L_i(x_j) should equal delta_ij');
        end
        
        function testInterpolationProperty(testCase)
            % Lagrange basis can exactly interpolate polynomial data
            % If we have values y at knots, interpolant is sum_i y_i * L_i(x)
            
            knots = testCase.basis.knots;
            degree = testCase.basis.degree;
            
            % Create test data: evaluates a polynomial at nodes
            test_poly_coeffs = rand(degree + 1, 1);
            
            % Evaluate polynomial at knots
            y_at_knots = polyval(flip(test_poly_coeffs), knots);
            
            % Evaluate Lagrange basis at a test point
            test_pts = linspace(testCase.basis.domain(1), testCase.basis.domain(2), 10);
            v = testCase.basis.evaluate(test_pts);
            
            % Interpolant using Lagrange basis
            interpolant = v * y_at_knots;
            
            % Expected: polynomial evaluated at test points
            expected = polyval(flip(test_poly_coeffs), test_pts)';
            
            testCase.verifyEqual(interpolant, expected, 'RelTol', 1e-12, ...
                'Lagrange interpolation property failed');
        end
        
        function testFirstDerivative(testCase)
            % Test that derivatives agree with finite differences

            test_pts = linspace(testCase.basis.domain(1), testCase.basis.domain(2), 8);
            delta = 1e-8;
            plus_pts=test_pts+delta;
            minus_pts=test_pts-delta;
            if ~testCase.basis.periodic
                plus_pts=min(plus_pts, testCase.basis.domain(2));
                minus_pts=max(minus_pts, testCase.basis.domain(1));
            end
            delta=plus_pts-minus_pts;

            % Evaluate function and first derivative
            [v1] = testCase.basis.evaluate(test_pts, parDers(1));
            % Finite difference approximation
            [v_plus] = testCase.basis.evaluate(plus_pts);
            [v_minus] = testCase.basis.evaluate(minus_pts);
            fd_deriv = (v_plus - v_minus) ./ delta(:);
            
            testCase.verifyEqual(v1, fd_deriv, 'RelTol', 1e-5, ...
                'Derivative values do not match finite difference approximation');
        end
        
        
        function testPartitionOfUnity(testCase)
            % Sum of all Lagrange basis functions should equal 1
            % (polynomial space property)
            
            test_pts = linspace(testCase.basis.domain(1), testCase.basis.domain(2), 20);
            [v] = testCase.basis.evaluate(test_pts);
            
            % Sum along basis dimension
            sums = sum(v, 2);
            
            % Should be close to 1
            expected = ones(npts(test_pts), 1);
            testCase.verifyEqual(sums, expected, 'RelTol', 1e-13, ...
                'Partition of unity property violated');
        end
        
        function testDimensionMatchesNodeCount(testCase)
            % Dimension should equal number of nodes
            num_nodes = numel(testCase.basis.knots);
            testCase.verifyEqual(testCase.basis.dimension, num_nodes, ...
                'Dimension should equal number of nodes');
            
            % Degree should be num_nodes - 1
            expected_degree = num_nodes - 1;
            testCase.verifyEqual(testCase.basis.degree, expected_degree, ...
                'Degree should be num_nodes - 1');
        end
        
        function testUniformNodeSpacing(testCase)
            % Test with uniformly spaced nodes (standard case)
            
            num_nodes = 4;
            knots = linspace(-1, 1, num_nodes);
            basis = lagrange(knots, [-1, 1], false);
            
            testCase.verifyEqual(basis.dimension, num_nodes);
            
            % Evaluate at nodes
            [v, a] = basis.evaluate(knots);
            v_full = a.toFull(v);
            
            % Should be identity
            testCase.verifyEqual(v_full, eye(num_nodes), 'AbsTol', 1e-13);
        end
        
        function testNonuniformNodeSpacing(testCase)
            % Test with non-uniformly spaced nodes
            
            % Chebyshev-like spacing
            num_nodes = 5;
            knots = cos(pi * (0:num_nodes-1) / (num_nodes - 1));
            basis = lagrange(knots, [-1, 1], false);
            
            % Should still satisfy Kronecker delta property
            [v, a] = basis.evaluate(knots);
            v_full = a.toFull(v);
            
            testCase.verifyEqual(v_full, eye(num_nodes), 'AbsTol', 1e-12, ...
                'Kronecker delta must hold for non-uniform nodes');
        end
        
    end
    
end
