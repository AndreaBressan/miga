function ncfs=oslo(ks, nks)
% oslo
% computes the coefficients the expansion of a coarse B-splines in finer knotvector
% ks are the local coarse knots for each fine knot in nks 
deg=size(ks,1)/2;
l=deg;
r=deg+1;
ncfs=zeros(deg+1,size(ks,2));
ncfs(1,:)=1;
for s=1:deg
    ncfs(s+1,:)=ncfs(s,:).*(nks(s+1,:)-ks(l,:))./(ks(l+s,:)-ks(l,:));
    for j=s-1:-1:1
        ncfs(j+1,:)=ncfs(j,:).*(nks(s+1,:)-ks(l+j-s,:))./(ks(l+j,:)-ks(l+j-s,:))...
            +ncfs(j+1,:).*(ks(r+j,:)-nks(s+1,:))./(ks(r+j,:)-ks(r+j-s,:));
    end
    ncfs(1,:)=ncfs(1,:).*(ks(r,:)-nks(s+1,:))./(ks(r,:)-ks(r-s,:));
end
end
