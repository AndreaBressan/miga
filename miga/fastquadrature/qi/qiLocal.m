function qi=qiLocal(basis,nodes)
% This is not local!!!
arguments
    basis (1,1) aBasis
    nodes (1,:) double
end
qi.basis=basis;
qi.nodes=nodes;
qi.mat=sparse(basis.dimension,npts(nodes));

[v,eval_sp]=basis.evaluate(nodes);
V=eval_sp.toSparse(v);

[pt, fun]=eval_sp.pos();

for b=1:basis.dimension
    lpt= fun==b & v~=0; 
    funPts = unique(pt(lpt));
    lfn=any(pt==funPts',2) & v~=0;
    funTst = unique(fun(lfn));
    rhs=zeros(numel(funTst),1);
    rhs(funTst==b)=1;
    mat=full(V(funPts,funTst)');
    b_op=lsqminnorm( mat , rhs );

    qi.mat(b,funPts)=b_op;
end
end