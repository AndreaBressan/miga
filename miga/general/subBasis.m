classdef subBasis < aBasis
    % subBasis is a linear combination in the space spanned by a set of
    % functions: a single function, a basis or generating set.

    properties(GetAccess=public,SetAccess=protected)
        basis (1,1)   aBasis=bSpline(0,[-inf,inf])
        coefs (:,:,:) double=double.empty
    end

    methods
        function obj = subBasis(basis,coefs,rangeSize)
            % coefs is a tensor of size
            % [basis.dimension,subBasis.dimension,prod(rangeSize)]
            arguments
                basis      (1,1)    aBasis
                coefs      (:,:,:)  double
                rangeSize  (:,1)    double = double.empty
            end
            assert( size(coefs,1)==basis.dimension, ...
                "The number of rows of coefficients must be equal to the dimension of the basis.")
            assert(size(coefs,3)==1||prod(basis.rangeSize)==1,...
                "Cannot compute vector-valued linear combination of vector-valued functions.")
            assert(isempty(rangeSize)||prod(rangeSize)==size(coefs,3)||prod(rangeSize)==prod(basis.rangeSize),...
                "Cannot compose a function of the required shape with the given basis and coefficients.")

            if isa(basis,"subBasis")
                obj.basis = basis.basis;
                obj.coefs = pagemtimes(basis.coefs,coefs);
            else
                obj.basis = basis;
                obj.coefs = coefs;
            end

            obj.domain = basis.domain;
            obj.periodic = basis.periodic;
            obj.dimension = size(coefs,2);
            if isempty(rangeSize)
                if size(coefs,3)>1
                    obj.rangeSize=size(coefs,3);
                elseif prod(basis.rangeSize)>1
                    obj.rangeSize=basis.rangeSize;
                else
                    obj.rangeSize=double.empty;
                end
            else
                obj.rangeSize=rangeSize;
            end
        end

        function [b]  = supports(obj)
            b = repmat(obj.domain,1,1,obj.dimension);
        end

        function [v,a] = evaluateCombinationI(obj,coefs,pts,ders)
            nC=tApplyOp(obj.coefs,coefs,2);
            [vt,at]=obj.basis.evaluateCombinationI(nC,pts,ders);
            [v,a]=obj.reshapeOutput(vt,at);
        end

        function [v,a] = reshapeOutput(obj,vt,at)
            rangeDim=prod(obj.rangeSize);
            basisRangeDim=prod(obj.basis.rangeSize);
            if rangeDim==basisRangeDim
                v=vt; a=at;
            else
                nPts=at.shape(1);
                if isa(at,'rSparsity')
                    [v,att]=obj.reshapeOutput(vt,at.original);
                    a=rSparsity(att,at.newRow,at.newShape);
                    return
                elseif isa(at,'gSparsity')
                    a=gSparsity(at.model*repelem(speye(obj.dimension),rangeDim,1)~=0);
                    [r,ct]=at.pos();
                    g=mod(ct-1,rangeSize)+1;
                    c=floor(ct-1/rangeSize)+1;
                    p=(a.nzPos(r,c)-1)*rangeDim+g;
                    reshapeF=@(w) modifyOutput(w,rangeDim,p,a.nnz);
                elseif isa(at,'fSparsity') % full
                    a=fSparsity([nPts,obj.dimension]);
                    shw=[basisRangeDim,nPts,obj.dimension,rangeDim/basisRangeDim];
                    shv=[rangeDim,1,nPts*obj.dimension];
                    reshapeF=@(w) reshape(permute(reshape(w,shw),[1,4,2,3]),shv);
                else
                    error('not implemented')
                end
                v=cellfun(reshapeF,vt,'UniformOutput',false);

            end
        end
    end

    methods
        function [v,a] = evaluateI(obj,pts,ders)
            [vt,at]=obj.basis.evaluateCombinationI(obj.coefs(:,:),pts,ders);
            [v,a]=obj.reshapeOutput(vt,at);
        end
    end

end


function nv=modifyOutput(v,rDim,pos,nDim)
    nv=zeros(rDim,1,nDim);
    nv(pos)=v;
end
