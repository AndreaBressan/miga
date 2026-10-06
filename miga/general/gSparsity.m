classdef gSparsity < aSparsity
    properties
        model
    end

    methods
        function obj=gSparsity(model)
            arguments
                model {mustBeSparse}
            end
            obj.nnz=nnz(model);
            obj.shape=size(model);
            [rs,cs,ord]=find(model);
            if islogical(model) || all(ord==1)
                obj.model=sparse(rs,cs,1:obj.nnz,obj.shape(1),obj.shape(2));
            else
                obj.model=model;
                assert(isequal(sort(ord(:)),(1:obj.nnz)'), ...
                    ['To construct a gSparsity use a logical sparse matrix'...
                     ' or a matrix where the entries specify the positions'...
                     ' in the vector of non-zeros']);
            end
        end
        
        function [out] = nnzColImpl(obj,col)
            out=full(sum(obj.model(:,col)~=0,1))';
        end

        function [out] = nnzColAllImpl(obj)
            out=full(sum(obj.model~=0,1))';
        end
        
        function [out] = nnzRowImpl(obj,row)
            out=full(sum(obj.model(row,:)~=0,2));
        end
        
        function [out] = nnzRowAllImpl(obj)
            out=full(sum(obj.model~=0,2));
        end

        function [nz] = nzPosImpl(obj,row,col)
            nz=full(obj.model(sub2ind(obj.shape,row,col)));
        end

        function [row,col] = posImpl(obj,nz)
            [row,col]=obj.posAllImpl();
            row=row(nz);
            col=col(nz);
        end

        function [row,col] = posAllImpl(obj)
            [row,col]=find(obj.model);
        end

        function out=toSparseImpl(obj,data,dir)
            sz=size(data,1:max(ndims(data),dir));
            [row,col,ord]=find(obj.model);
            rep=prod(sz(1:(dir-1)));
            if rep>1
                row=(1:rep)'+rep*(row(:)'-1);
                col=repelem(col(:),rep);
                ord=(1:rep)'+rep*(ord(:)'-1);
            end
            out=sparse(row(:),col,reshape(data(ord),[],1),obj.shape(1)*rep,obj.shape(2));
        end
    end
end