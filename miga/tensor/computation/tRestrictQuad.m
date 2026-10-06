function Qside=tRestrictQuad(Q,side)
arguments
    Q    (1,1) tQuadrature
    side (:,1) double 
end 
Qside = tQuadrature(Q.components(side==0));
end