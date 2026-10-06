function [QT,R,P,D] = pageQR(J,part)
%PAGEQR decompose an array of matrices using Q R decomposition
% On each page k the results are such that
%
% J(:,:,k)=QT(:,:,k)'*R(:,P,k)
%
% The additional output D contains the volume transformation of J for the
% specified part.
%
% The implementation is tailored to invert the jacobian J of a geometric
% map G defined on a box. The optional argument 'part' specifies which part
% of the box should be considered the domain of G in such a way that the
% normal are available after the QR decomposition.

% A zero in part means that the corresponding direction is free, a one that
% the corresponding variable is fixed to the maximum of the box interval,
% a -1 means the minimum of the box interval, e.g. [0,1] is the top segment
% of the boundary of a 2D box, while [-1,1] is the top-left corner.

% The first group of rows of QT, corresponding zeros entries in part,
% and are an orthonormal basis of the tangent space to the image by G of 
% the specified part. They are set to have a positive scalar product with
% the corresponding partial derivatives of G, i.e. R has positive diagonal
% entries.
%
% The second group of rows corresponds to the non-zeros entries of part
% and complements the first group to an orthonormal basis of the tangent
% to the image by G of the whole box. The entries have positive or negative
% scalar product with the i-th column of J respectively if part is -1 or 1,
% so that they represent a basis of the outer-normals.
% For the application of the Stokes' theorem on the integration of
% differential forms it is thus necessary to keep track of the sign due to
% order of the dx_i outside of this function.
%
% The third group complements the basis to the whole space, it non empty
% then the determinant of QT is positive.

% TODO add pivoting capabilities
rs=size(J,1);
cs=size(J,2);
nm=size(J,3:max(3,ndims(J)));
num=prod(nm);
assert(rs>=cs);

% init permutation of directions
if nargin==1; part=zeros(1,cs); end
P=[find(part==0);find(part~=0)];
nfree=sum(part==0);
Ppart=part(P);
% init data
QT=repmat(eye(rs),[1,1,nm]);
R=reshape(J(:,P,:),size(J));
D=ones([1,1,prod(nm)]);
if rs==1
    neg=R<0;
    if nfree==1; D=R; else, D=ones(size(R)); end
    R=abs(R);
    QT=reshape(-neg+~neg, size(R));
    return
end

for c=1:min(cs,rs-1)
    v= R(c:end,c,:);
    s= (v(1,:,:)>=0)-(v(1,:,:)<0);
    n=vecnorm(v,2,1);
    v(1,:,:)=v(1,:,:)+n.*s;
    v=v./vecnorm(v,2,1);
    % compute update
    RU=pagemtimes(v,pagemtimes(2*v,'transpose',R(c:end,c:end,:),'none'));
    QU=pagemtimes(v,pagemtimes(2*v,'transpose',QT(c:end,:,:),'none'));
    % apply Hausholder reflection
    R(c:end,c:end,:)=reshape(-s,[1,1,num]).*(R(c:end,c:end,:)-RU(:,:,:));
    QT(c:end,:,:)   =reshape(-s,[1,1,num]).*(QT(c:end,:,:)-QU(:,:,:));
    % stop propagation of truncation errors
    R(c+1:end,c,:)=0; 
    % update determinants
    if c<=nfree;  D=D.*R(c,c,:); end
    if nfree==rs; D= -D.*(-s).^mod(rs-c+1,2);end
end
if c<cs % J is square, c+1=rs=cs
    neg=R(cs,cs,:)<0;
    s= -neg +~neg;
    if cs>nfree; s=s*Ppart(cs); end
    R(cs,:,:)=s.*R(cs,:,:);
    QT(cs,:,:)=s.*QT(cs,:,:);
    if cs==nfree; D=D.*R(cs,cs,:); end  
    if nfree==rs; D=D.*s;end 
end
end

