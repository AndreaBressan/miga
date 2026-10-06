%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example of solution od advection-diffusion-reaction equation on a       %
% single-patch geometry.                                                  %
%                                                                         %
%  -div(K*grad(u) + grad(u)*[adv1;adv2] + react*u = f in geo([0,1]^n)     %
%  u         = uD                                           in GammaD     %
%  d_n u     = uN                                           in GammaN     %
%  u + d_n u = uR                                           in Omega      %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Problem description
problemData.geo = quarterRing(1,2);
problemData.dim = size(problemData.geo.domain,1);
domain = repmat([-inf,inf],[problemData.dim,1]);
problemData.drchltSides = [-1 0]';
problemData.nmnnSides   = [1 0]';
problemData.rbnSides    = [0 -1;0 1]';
syms 'x1' 'x2';
K11 = 1;
K12 = 0;
K21 = K12;
K22 = 1;
adv1 = 1;
adv2 = 1;
react = 1;

% manufactured solution
uEx = -sin(2*atan(x2/x1))*sin(pi*sqrt((x1^2+x2^2)));
uLapEx = diff(K11*diff(uEx,x1)+K12*diff(uEx,x2),x1)+diff(K21*diff(uEx,x1)+K22*diff(uEx,x2),x2);
uAdv = adv1*diff(uEx,x1,1)+adv2*diff(uEx,x2,1);
uReact = react*uEx;
uNEx = (x1*diff(uEx,x1)+x2*diff(uEx,x2))/sqrt(x1^2+x2^2);% fitted for the current geometry (QoR)
uREx1 = simplify(-diff(uEx,x2)+uEx);% fitted for the current geometry (QoR)
uREx2 = simplify(-diff(uEx,x1)+uEx);% fitted for the current geometry (QoR)
% make functions out of expressions
problemData.eqCoefs = exprBasis(domain,[string(K11),string(K12),string(K21),string(K22),string(adv1),string(adv2),string(react)]);
uSol = exprBasis(domain,string(uEx));
problemData.source = exprBasis(domain,string(-uLapEx+uAdv+uReact));
problemData.uD = uSol;
problemData.uN = exprBasis(domain,string(uNEx));
problemData.uR(1) = exprBasis(domain,string(uREx1));
problemData.uR(2) = exprBasis(domain,string(uREx2));

%% Discretization Parameters
methodData.deg  = [ones(1,problemData.dim)]*3;
methodData.nsub = [ones(1,problemData.dim)]*20;

%% Solver call
sol = solveADReq(problemData,methodData);

%% Post-processing
uh_val = sol.uh.evaluate(sol.Q.nodes);
u_val = uSol.evaluate(sol.geo.evaluate(sol.Q.nodes));
err_val=uh_val-u_val;
errL2 = sqrt(sol.Q.integrate(H1pullback(sol.geo,preEvaluatedFunction(err_val.^2,domain))));

%% Plots
phyPts = sol.geo.evaluate(sol.Q.nodes);
subplot(1,3,1); 
    surfLines(reshape(phyPts(1,:),cellfun(@numel,sol.Q.nodes)),reshape(phyPts(2,:),cellfun(@numel,sol.Q.nodes)),reshape(uh_val,cellfun(@numel,sol.Q.nodes)));
    title('uh');xlabel('x'); ylabel('y');
    view(3); axis equal;
subplot(1,3,2); 
    surfLines(reshape(phyPts(1,:),cellfun(@numel,sol.Q.nodes)),reshape(phyPts(2,:),cellfun(@numel,sol.Q.nodes)),reshape(u_val,cellfun(@numel,sol.Q.nodes)));
    title('uSol');xlabel('x'); ylabel('y');
    view(3); axis equal;
subplot(1,3,3);
    plotGeometry(sol.geo);
    surfLines(reshape(phyPts(1,:),cellfun(@numel,sol.Q.nodes)),reshape(phyPts(2,:),cellfun(@numel,sol.Q.nodes)),reshape(err_val,cellfun(@numel,sol.Q.nodes)));
    daspect ([1 1 max(abs(err_val),[],'all')]);
    title('Err');xlabel('x'); ylabel('y');  
    view(3); 
