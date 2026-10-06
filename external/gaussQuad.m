function [m] = gaussQuad(n, name, alpha, beta)
% namedQuads wraps Chebfun computation of n quadrature nodes and weights
% for the following quadratures
% .Legendre, the standard Gauss quadrature;
% .Jacobi, the Gauss quadrature for the weight (x+1)^alpha*(x-1)^beta;
% .Lobatto, modified Gauss-Jacobi quadrature that includes the endpoints;
% .Radau, modified Gauss quadrature that includes the first endpoint;
% .Chebychev, a special case of Gauss-Jacobi quadrature;
% .Clenshaw-Curtis, a modification of Chebychev quadrature for adaptive
%            integration where some knots are reused after subdivision
arguments
    n     (1,1) double {isinteger} = 0
    name  (1,1) string = 'legendre'
    alpha (1,1) double = 0
    beta  (1,1) double = 0
end
switch lower(name)
    case 'legendre'
        [nodes,weights] = legpts(n);
        % weight is 1 on [-1,1]
        % exact up to degree 2*n-1
%    case 'hermite'
%        [nodes,weights] = hermpts(n);
%        % weight is e^(-x^2) on R
%        % exact up to degree 2*n-1
%    case 'laguerre'
%        [nodes,weights] = lagpts(n);
%        % weight is e^(-x) on R>0
%        % exact up to degree 2*n-1
    case 'jacobi'
        [nodes,weights] = jacpts(n,alpha,beta);
        % weight is (x+1)^alpha*(x-1)^beta on [-1,1]
        % exact up to degree 2*n-1
    case 'lobatto'
        [nodes,weights] = lobpts(n, alpha, beta);
        % weight is (x+1)^alpha*(x-1)^beta on [-1,1]
        % but the endpoints are quadrature nodes
        % exact up to degree 2*n-3
    case 'chebychev'
        [nodes,weights] = jacpts(N, -.5, -.5);
    case 'clenshaw-curtis'
        [nodes,weights] = chebpts(n);
    case 'radau'
        [nodes,weights] = radaupts(n);
        % includes -1 and is exact up to degree 2*n-3
    otherwise
        error('unknown qudrature name')
end
m=struct("name",lower(name),"weights",weights(:)',"nodes",nodes(:)');
end

