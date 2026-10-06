function out=tApplyOp(t,op,dir)
% tApplyOp
% applies the operator op to the nD array t in the direction dir
arguments
    t  double
    op
    dir (1,1) double = 1
end
isz=size(t,1:max(dir,ndims(t)));
osz=isz;
osz(dir)=size(op,2);
sb=prod(isz(1:dir-1));
sa=prod(isz(dir+1:end));

if ~issparse(op) || ( (numel(op)/nnz(op)<30) && numel(op)<2^24)
    op=full(op);
    shpt=reshape(t,[prod(isz(1:dir-1)),isz(dir),isz(dir+1:end)]);
    tmp=pagemtimes(shpt,op);
elseif sb==1
    tmp=op'*reshape(t,isz(dir),[]);
elseif sa==1
    tmp=reshape(t,[],isz(dir))*op;
else
    tmp=zeros(sb,osz(dir),sa);
    shpt=reshape(t,sb,isz(dir),sa);
    for it=1:sa
        tmp(:,:,it)=shpt(:,:,it)*op;
    end
end
out=reshape(tmp,osz);
end