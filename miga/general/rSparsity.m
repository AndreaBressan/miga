classdef rSparsity < aSparsity
    properties
        original (1,1) aSparsity = fSparsity([0,0])
        newRows  (:,1) double
        oldRows  (:,1) double
    end

    methods
        function obj=rSparsity(original,newRow,newShape)
            arguments
                original (1,1) aSparsity
                newRow   (:,1) double
                newShape (1,2) double
            end
            obj.original=original;
            obj.nnz=original.nnz;
            obj.shape=newShape;
            obj.newRows=newRow;
            obj.oldRows(newRow)=1:original.shape(1);
            assert(all(newRow<=newShape(1)),'wrong row positions')

        end
        
        function [out] = nnzColImpl(obj,col)
            out=obj.original.nnzColImpl(col);
        end

        function [out] = nnzColAllImpl(obj)
            out=obj.original.nnzColAllImpl();
        end
        
        function [out] = nnzRowImpl(obj,row)
            oldR=obj.oldRows(row);
            out=zeros(numel(row),1);
            out(oldR~=0)=obj.original.nnzRowImpl(oldR~=0);
        end
        
        function [out] = nnzRowAllImpl(obj)
            out(obj.oldRows~=0)=obj.original.nnzRowAllImpl();
        end

        function [nz] = nzPosImpl(obj,row,col)
            nz=obj.original.nzPosImpl(obj.oldRows(row),col);
        end

        function [row,col] = posImpl(obj,nz)
            [rowT,col] = obj.original.posImpl(nz);
            row=obj.newRows(rowT);
        end

        function [row,col] = posAllImpl(obj)
            [rowT,col] = obj.original.posAllImpl();
            row=obj.newRows(rowT);       
        end
    end
end