# miga
A matlab implementation of isogeometric analysis based on sum factorization

## Installation

1. Clone the repository and initialize the submodules[^1] by using the shell command
```bash
 git clone --recurse-submodules git@github.com:AndreaBressan/miga.git REPOSITORY_PATH
```
2. Start matlab and add the repository to its path: in the command window type
```matlab
 addpath(genpath(REPOSITORY_PATH));
```
3. Optionally compile the provided mesh file to speed-up sparse matrix forming[^2]: in matlab command window type
```matlab
 cd REPOSITORY_PATH
 mex_build_all
```

[^1]: `miga` uses [chebfun](https://github.com/chebfun/chebfun) to compute quadrature weights and nodes.
[^2]: The non-zeros computed by the sum-factorization algorithm are positioned in a known order so that sorting them requires order ($O(\textrm{nnz})$) operations instead of order ($O(\textrm{nnz} \log(\textrm{nnz}))$) required by calling `sparse(row,col,val,n,m)`.

## Library structure

The library is built around three interfaces: basis, integration, and sparsity.

1. Bases represent functions or sets of functions, they can be evaluated on points and on grids. In the following example basis is a set of function (the B-spline basis), spline is a set of a single function (a spline curve).
```matlab
degree=3;
knots=[1 1 1 1 2 3 4 5];
basis=bSpline(degree,knots);
xs=linspace(1,5,100);
ys=basis.evaluate(xs); % evaluate values
ds=basis.evaluate(xs,parDers(1)); % evaluate first derivative
basis.plot(); % plot all basis functions
coefficients=rand(basis.dimension,1);
f=subBasis(basis,coefficients);
f.plot();
```
2. Integration represents a quadrature rule on the whole domain. It can integrate functions, it can integrate products of a function times the basis functions or a bilinear form with test and trial spaces (i.e. a Petrov-Galerkin matrix)
```
elements=bases2elements(basis);
quadrature=uElementQuadrature(elements, gaussQuad(degree+1));
fIntegral=quadrature.integrate(f); % integrate the provided function
test=basis;
trial=basis;
fTest=quadrature.integrate(f,test); % integrate all products f basis_function
mass=quadrature.integrate(constantFunction(1,test.domain),test,zeros(1,1,2),trial); % integrates all pairs of test and trial with constant weight 1
% replacing zeros(1,1,2) with reshape([1,2],1,1,2) would integrate the first derivative of the test against the second derivative of the trial.
max(abs(coefficients-mass\fTest)) % print the error of the solver of the L^2 projection
```
3. Sparsity represents the structure on non-zeros in the Galerkin matrices.

## Poisson with homogeneous Dirichlet conditions
```matlab
% domain
geo= quarterRing(1,2); 
domain=geo.domain;
dim=size(domain,1);

% right hand side
R2=repmat([-inf,inf],2,1); % the domain is the whole plane
physicalF=exprBasis(R2,"sin(3*x1+7*x2^2)"); % a function defined on the plane with a given expression

% describe the operator as pairs of derivatives and coefficients
% column (:,i,1) is the derivative of the test in the i-th direction
% and is integrated multiplied by the derivative of the trial in the 
% direction column(:,i,2)
% the coefficient function for each pair of derivatives is passed as
% a basis type, in this case a set of dim constant functions all equal to 1. 
physicalPairs=repmat(eye(dim),1,1,2); 
physicalCoefs=constantFunction(ones(dim,1),R2,[],dim); 

% set test and trial space splines that have zero value at the boundary
degree=repmat(3,1,dim);;
n=repmat(4,1,dim);

for i=1:dim
	knots=[domain(i,1)*ones(1,degree(i)-1),linspace(domain(i,1),domain(i,2),n(i)+1),domain(i,2)*ones(1,degree(i)-1)]; 
	uBasis(i)=bSpline(degree(i),knots); % one can have different degrees and knots per direction
end
basis=tBasis(uBasis);
test=basis;
trial=basis; % use a different basis for Petrov-Galerkin cases

% choose integration
uQuad=arrayfun(@(i) uElementQuadrature(bases2elements(uBasis(i)), gaussQuad(degree(i)+1)),1:dim);
quad=tQuadrature(uQuad);

% map operator and rhs to the parametric domain
parametricCoefs=H1pullback(geo,physicalCoefs,physicalPairs);
parametricF=H1pullback(geo,physicalF);

% integrate
A=quad.integrate(parametricCoefs,test,parametricCoefs.newRequest,trial);
rhs=quad.integrate(parametricF,test);
uCoefs=A\rhs;
u=subBasis(trial,uCoefs);

% evaluate and plot the solution
NumPoints=100*ones(1,dim);
parametricXs=arrayfun(@(i)linspace(domain(i,1),domain(i,2),NumPoints(i)),1:dim,'UniformOutput',false);
ys=u.evaluate(parametricXs);
physicalXs=geo.evaluate(parametricXs);

surfLines(reshape(physicalXs(1,:),NumPoints),reshape(physicalXs(2,:),NumPoints),reshape(ys,NumPoints));
```

