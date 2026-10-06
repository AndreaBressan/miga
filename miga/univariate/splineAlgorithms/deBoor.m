function v=deBoor(der,x,lks)
nPts=size(lks(:,:),2);
deg=size(lks,1)/2;
v=zeros(deg+1,nPts,size(lks,3));
if der>deg
    return
end
st=1; for s=1:der; st=st*(deg+1-s); end
v(1,:)=st;
l=deg;
for s=1:deg-der
    J=1:s-1;
    v(s+1,:,:)=v(s,:,:).*(x-lks(l,:,:))./(lks(l+s,:,:)-lks(l,:,:));
    v(J+1,:,:)=v(J,:,:).*(x-lks(l+J-s,:,:))./(lks(l+J,:,:)-lks(l+J-s,:,:))...
               +v(J+1,:,:).*(lks(l+1+J,:,:)-x)./(lks(l+1+J,:,:)-lks(l+1+J-s,:,:));
    v(  1,:,:)=v(1,:,:).*(lks(l+1,:,:)-x)./(lks(l+1,:,:)-lks(l+1-s,:,:));
end
for s=deg-der+1:deg
    J=1:s-1;
    v(s+1,:,:)= v(s,:,:)./(lks(l+s,:,:)-lks(l,:,:));
    v(J+1,:,:)= v(J,:,:)./(lks(l+J,:,:)-lks(l+J-s,:,:))-v(J+1,:,:)./(lks(l+1+J,:,:)-lks(l+1+J-s,:,:));
    v(  1,:,:)=-v(1,:,:)./(lks(l+1,:,:)-lks(l+1-s,:,:));
end
end