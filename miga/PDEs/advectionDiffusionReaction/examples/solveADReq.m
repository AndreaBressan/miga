function sol = solveADReq(problemData,methodData)

%% Unwrapping input structures
uN          =problemData.uN;
uD          =problemData.uD;
uR          =problemData.uR;
geo         =problemData.geo;
drchltSides =problemData.drchltSides;
nmnnSides   =problemData.nmnnSides;
rbnSides    =problemData.rbnSides;
eqCoefs     =problemData.eqCoefs;
source      =problemData.source;

deg         =methodData.deg;
nsub        =methodData.nsub;

dim = size(geo.domain,1);

%% Assembling
for iDim = dim:-1:1
    knots = [zeros(1,deg(iDim)),linspace(0,1,nsub(iDim)+1),ones(1,deg(iDim))];
    uB(iDim) = bSpline(deg(iDim),knots);
end
B = tBasis(uB);
uEls = arrayfun(@(uBas)bases2elements(uBas),uB,'UniformOutput',false);
uQ = arrayfun(@(iDim)uElementQuadrature(uEls{iDim},gaussQuad(deg(iDim)+1)),1:dim);
Q = tQuadrature(uQ);
gradgrad = cat(3,repmat(eye(dim),[1,dim]),repelem(eye(dim),1,dim));
valAdv = cat(3,zeros(dim),eye(2));
valval = cat(3,zeros(dim,1),zeros(dim,1));
pullback = H1pullback(geo,eqCoefs,cat(2,gradgrad,valAdv,valval));
Mat = Q.integrate(pullback,B,pullback.newRequest,B);
Rhs = Q.integrate(H1pullback(geo,source),B);
u = zeros(B.dimension,1);

%% Imposition of Neumann boundary conditions
for iSide = 1:size(nmnnSides,2)
    Rhs = Rhs+tSideQuadrature(Q,geo.domain,nmnnSides(:,iSide)).integrate(H1pullback(geo,uN,[],nmnnSides(:,iSide)),B);
end

%% Imposition of Robin boundary conditions
for iSide = 1:size(rbnSides,2)
    Qside = tSideQuadrature(Q,geo.domain,rbnSides(:,iSide));
    valval = H1pullback(geo,constantFunction(1,repmat([-inf,inf],dim,1)),[],rbnSides(:,iSide));
    Mat = Mat + Qside.integrate(valval,B,[],B);
    Rhs = Rhs + Qside.integrate(H1pullback(geo,uR(iSide),[],rbnSides(:,iSide)),B);
end

%% Projection the Dirichlet datum for each side
% compute the average for dofs appearing in more than one side
% (sideCount is the count of the number of sides in which every dof appears)
sideCount = zeros(size(u));
for iSide = 1:size(drchltSides,2)
    [~,idx] = tRestrict(B,drchltSides(:,iSide));
    Qside = tSideQuadrature(Q,geo.domain,drchltSides(:,iSide));
    valval = H1pullback(geo,constantFunction(1,repmat([-inf,inf],[dim,1])),[],drchltSides(:,iSide));
    Mside = Qside.integrate(valval,B,[],B);
    DrSide = Qside.integrate(H1pullback(geo,uD,[],drchltSides(:,iSide)),B);
    u(idx) = u(idx)+Mside(idx,idx)\DrSide(idx);
    sideCount(idx) = sideCount(idx)+1;
end
intDofs = sideCount==0;
sideCount = max(sideCount,1);
u = u./sideCount;
Rhs = Rhs-Mat*u;

%% Linear solver
u(intDofs) = Mat(intDofs,intDofs)\Rhs(intDofs);

%% Preparing output
sol.uh = subBasis(B,u,[]);
sol.geo = geo;
sol.Q = Q;

end