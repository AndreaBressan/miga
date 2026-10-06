function qi=qiGreville(basis)
% THis is not local!!!
arguments
    basis (1,1) bSpline
end
qi.basis=basis;
qi.nodes=basis.greville;
qi.mat=inv(qi.basis.evaluate(qi.nodes(:)'));
end