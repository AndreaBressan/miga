
function sol=solveLS(deriv,coefs,Q,basis,options)
arguments
    deriv                   (:,:,1) double       {mustBeNonempty}
    coefs                   (1,1)   aBasis       {validateCoefs(coefs,deriv)}
    Q                       (1,1)   aIntegration
    basis                   (1,1)   aBasis
    options.rhs             (:,1)   aBasis       = aBasis.empty
    options.dirichlet       (:,1)   aBasis       = aBasis.empty
    options.dirichlet_sides (:,:)   double       = double.empty 
    options.neumann         (:,1)   aBasis       = aBasis.empty
    options.neumann_sides   (:,:)   double       = double.empty 
    options.robin           (:,1)   aBasis       = aBasis.empty
    options.robin_sides     (:,:)   double       = double.empty
    options.geometry        (:,1)   aBasis       = aBasis.empty
end
domain=basis.domain;


dim=size(Q.domain,1);
allBases=[coefs;basis;options.dirichlet;options.neumann;options.robin];
allDim=arrayfun(@(b)size(b.domain,1),allBases);
assert(all(allDim==dim), "Incompatible dimensions, check the dimensions of the domains of coefficients, basis, quadrature, boundary conditions");

LS_deriv = cat(3,repelem(deriv,1,size(deriv,2)),repmat(deriv,1,size(deriv,2)));
LS_coefs = productBasis(coefs,coefs);
RHS_deriv=deriv;
RHS_coefs=coefs;


if ~isempty(options.geometry)
    % TODO
    % handle Neumann and Robin conditions
    LS_coefs  = H1pullback(options.geometry,LS_coefs,LS_deriv);
    LS_deriv  = LS_coefs.newRequest;
    RHS_coefs = H1pullback(options.geometry,RHS_coefs,RHS_deriv);
    RHS_deriv = RHS_coefs.newRequest;
    if ~isempty(options.dirichlet)
        options.dirichlet=H1pullback(options.geometry,options.dirichlet);
    end
    if ~isempty([options.neumann_sides,options.robin_sides])
        error("Neumann and Robin conditions are not implemented when a geometry is present")
    end
end


% interior integrals
mat=Q.integrate(LS_coefs,basis,LS_deriv,basis);
mat=(mat+mat')/2; % very small errors can change the solver selection in matlab
if isempty(options.rhs)
    rhs=zeros(size(mat,1),1);
else
    rhs=Q.integrate(productBasis(RHS_coefs,options.rhs),basis,RHS_deriv);
end

% Dirichlet condition
u=zeros(size(mat,2),1);
sideCount=zeros(size(u));
for side = options.dirichlet_sides 
    [S_basis, idx] = tRestrict(basis, side);
    S_dirichlet    = tRestrict(options.dirichlet,side);
    S_quad         = tRestrictQuad(Q,side);
    S_mat          = S_quad.integrate(constantFunction(1,basis.domain(side==0,:)),S_basis,[],S_basis);
    S_rhs          = S_quad.integrate(S_dirichlet,S_basis);
    u(idx)=u(idx)+S_mat \ S_rhs;
    sideCount(idx)=sideCount(idx)+1;
end
intDofs   = sideCount==0;
sideCount = max(sideCount,1);
u=u./sideCount;
rhs=rhs-mat*u;

% Neumann condition
sideCount=zeros(size(u));
N_mat=sparse(basis.dimension,basis.dimension);
N_rhs=zeros(basis.dimension,1);
value   = zeros(dim,1);
for side = options.neumann_sides 
    normal  = side;
    N_quad  = tSideQuadrature(Q,domain,side);
    N_mat   = N_mat + N_quad.integrate(constantFunction(1,domain),basis,cat(3,value,normal),basis);
    N_rhs   = N_rhs + N_quad.integrate(options.neumann,basis);
    sideCount = sideCount+any(N_mat~=0,2);
end
neuman_idx = sideCount~=0;
neuman_count = max(sideCount(neuman_idx,1));
N_mat      = N_mat(neuman_idx,:)./neuman_count;
N_rhs      = N_rhs(neuman_idx)./neuman_count-N_mat*u;

sideCount=zeros(size(u));
R_mat=sparse(basis.dimension,basis.dimension);
R_rhs=zeros(basis.dimension,1);
for side = options.robin_sides 
    normal  = side;
    R_quad = tSideQuadrature(Q,domain,side);
    R_mat  = R_mat + R_quad.integrate(constantFunction([1,1],domain,[],2),basis,cat(3,[value,value],[normal,value]),basis);
    R_rhs  = R_rhs + R_quad.integrate(options.robin,basis);
    sideCount = sideCount+any(R_mat~=0,2);
end
robin_idx = sideCount~=0;
robin_count = max(sideCount(robin_idx,1));
R_mat     = R_mat(robin_idx,:)./robin_count;
R_rhs     = R_rhs(robin_idx)./robin_count-R_mat*u;

% restrict to inner dofs
mat=mat(intDofs,intDofs);
N_mat=N_mat(:,intDofs);
R_mat=R_mat(:,intDofs);
rhs=rhs(intDofs);

BC_mat=[N_mat;R_mat];
BC_num=size(BC_mat,1);

% total matrix
T_mat= [mat, BC_mat'; ...
        BC_mat, zeros(BC_num,BC_num)];
T_rhs=[ rhs;N_rhs;R_rhs];

% solve
%spparms('spumoni',1)
coefs=T_mat\T_rhs;
%spparms('spumoni',0)
u(intDofs)=coefs(1:size(mat,1));
sol=subBasis(basis,u,[]);
end

function validateCoefs(coefs,derivatives) 
if coefs.dimension~=size(derivatives,2)
    error('The number of derivative pairs and coeffients must match.') 
end 
end 