function [u, source, uD, uN, uR, coefs] = adrSolManufacturer(dim, sol, eqCoefs)

arguments
    dim    (1,1)  double
    sol    (1,1)  sym
    eqCoefs  (:,1)  sym % size(coefs,1) = dim^2 + dim + 1
end

coords = arrayfun(@(i) ['x', num2str(i)], 1:dim, 'UniformOutput', false);
syms(coords{:});
coords = sym(coords);

uGrad = sym(zeros(dim,1));
for iDim = 1:dim
    uGrad(iDim) = diff(sol,coords(iDim));
end
diffusion_uGrad = reshape(eqCoefs(1:dim^2),dim,dim)*uGrad;

uLapl = 0;
uAdv = 0;
for iDim = 1:dim
    uLapl = uLapl + diff(diffusion_uGrad(iDim),coords(iDim));
    uAdv = uAdv + eqCoefs(dim^2+iDim)*diff(sol,coords(iDim));
end

uReact = eqCoefs(end)*sol;

domain = repmat([-inf,inf],[dim,1]);

u = exprBasis(domain,string(sol));
source = exprBasis(domain,string(-uLapl+uAdv+uReact));
uD = u;
uN = []; % require geometry and nmnnSides
uR = []; % require geometry and rbnSides
coefs = exprBasis(domain,string(eqCoefs));

end