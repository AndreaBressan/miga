classdef testPoisson < aMigaTestCase

methods(Test)
function examplePoisson(testCase)

%% Problem description
dim = 2;
domain=repmat([0,1],[dim,1]);
drchltSides = [-1 0;1 0]';
nmnnSides   = [0,-1]';
rbnSides    = [0,1]';
thermalDiff = exprBasis(domain,repmat("1",1,dim)); % for later this must be a set of 2 functions
syms 'x1' 'x2';
% manufacture solution
uEx = (1+ x1.*x2+x1.^3 );%+exp(x1+x2.^2));
uLapEx =  diff(uEx,x1,2)+diff(uEx,x2,2);
uNEx   =  (2*x2-1).*diff(uEx,x2);     
uREx   =  (2*x2-1).*diff(uEx,x2)+uEx; 
% make functions out of expressions
uSol   = exprBasis(domain,string(uEx));
source = exprBasis(domain,string(-uLapEx));
uD     = uSol;
uN     = exprBasis(domain,string(uNEx));
uR     = exprBasis(domain,string(uREx));

%% Discretization parameters
deg = 3;
nsub = 10;

%% Assembling
knots = [zeros(1,deg),linspace(0,1,nsub+1),ones(1,deg)];
uB(1:dim) = bSpline(deg,knots);
B = tBasis(uB);
uEls = arrayfun(@(uBas)bases2elements(uBas),uB,'UniformOutput',false);
uQ = cellfun(@(uEls)uElementQuadrature(uEls,gaussQuad(deg+1)),uEls);
Q = tQuadrature(uQ);
gradgrad = repmat(eye(dim),[1,1,2]);
Mat = Q.integrate(thermalDiff,B,gradgrad,B);
Rhs = Q.integrate(source,B);
u = zeros(B.dimension,1);

for side = nmnnSides
    Rhs=Rhs+tSideQuadrature(Q,domain,side).integrate(uN,B);
end

for side = rbnSides
    Qside=tSideQuadrature(Q,domain,side);
    Mat  = Mat + Qside.integrate(constantFunction(1,domain),B,[],B);
    Rhs  = Rhs + Qside.integrate(uR,B);
end

% project the Dirichlet datum for each side and compute the average for dofs
% appearing in more than one side (sideCount is the count of the number of sides in
% which every dof appears)
sideCount=zeros(size(u));
for side = drchltSides 
    [Bs,idx] = tRestrict(B, side);
    [uDs]    = tRestrict(uD,side);
    Qside    = tRestrictQuad(Q,side);
    
    Mside=Qside.integrate(constantFunction(1,domain(side==0,:)),Bs,[],Bs);
    DrSide=Qside.integrate(uDs,Bs);
    u(idx)=u(idx)+Mside\DrSide;
    sideCount(idx)=sideCount(idx)+1;
end
intDofs=sideCount==0;
sideCount=max(sideCount,1);
u=u./sideCount;
Rhs=Rhs-Mat*u;

u(intDofs) = Mat(intDofs,intDofs)\Rhs(intDofs);
uf=subBasis(B,u,[]);

uh_val=uf.evaluate(Q.nodes);
u_val=uSol.evaluate(Q.nodes);
errL2=sqrt(Q.integrate(preEvaluatedFunction((uh_val-u_val).^2,domain)));
testCase.verifyEqual(errL2,0,"AbsTol",1e-10);

import matlab.unittest.diagnostics.Diagnostic;
import matlab.unittest.diagnostics.FigureDiagnostic;
f=testCase.invisibleFigure();


[X,Y]=ndgrid(Q.nodes{:});
subplot(1,3,1); uf.plot; title('computed solution'); xlabel('x'); ylabel('y');
subplot(1,3,2); uSol.plot; title('correct solution');xlabel('x'); ylabel('y');
subplot(1,3,3); surf(X,Y,reshape(uf.evaluate(Q.nodes)-uSol.evaluate(Q.nodes),size(X))); title('difference');xlabel('x'); ylabel('y');

testCase.log(3,FigureDiagnostic(f));
end
end
end

