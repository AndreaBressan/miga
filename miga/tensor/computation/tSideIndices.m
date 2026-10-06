function idx=tSideIndices(basis,side)
arguments
    basis (1,1) tBasis
    side  (:,1) double
end
dim=size(basis.domain,1);
sz=basis.vectorDimension;
free = (side==0);
rep=[max(sz.*free(:)',1),1];

ui=cell(dim,1);
for d=1:dim
    if free(d)
        rep(d)=1;
        ui{d}=repmat(reshape(1:sz(d),[ones(1,d-1),sz(d),1]),rep);
        rep(d)=sz(d);
    else
        ui{d}=repmat(sz(d)*(side(d)==1)+(side(d)==-1),rep);
    end
end
idx=sub2ind(sz,ui{:});
end