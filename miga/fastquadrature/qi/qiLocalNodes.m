function nodes=qiLocalNodes(basis)
arguments
    basis (1,1) bSpline
end
K=basis.knots;
start=[true;K(2:end)-K(1:end-1)>0];
startP=find(start);
uK=K(startP);
reps=diff(startP);reps=[reps; numel(K)-startP(end)+1];

nodesPerEl= ceil((reps(1:end-1)+reps(2:end))/2);
nodesPerEl(1)=basis.degree+1;
nodesPerEl(end)=basis.degree+1;
nodesPerEl=nodesPerEl+1;

elSize=uK(2:end)-uK(1:end-1);
elNodes=arrayfun(@(el) uK(el)+elSize(el)*(0:nodesPerEl(el)-1)/nodesPerEl(el),1:numel(elSize),"UniformOutput",false);

nodes=cat(2,elNodes{:},uK(end));

end