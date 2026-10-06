
classdef fSparsity <aSparsity

    methods
        function obj=fSparsity(shape)
            arguments
                shape (1,:) double
            end
            obj.nnz=prod(shape);
            obj.shape=shape;
        end
        function [out]=toFullImpl(obj,data,dir)
            out=reshape(data,[size(data,1:(dir-1)),obj.shape,size(data,(dir+1):ndims(data))]);
        end
        function [out] = nnzColImpl(obj,col)
            out=obj.shape(1)*ones(numel(col),1);
        end
        function [out] = nnzRowImpl(obj,row)
            out=obj.shape(2)*ones(numel(row),1);
        end
        function [row,col] = posImpl(obj,nz)
            [row,col]=ind2sub(obj.shape,nz);
        end
        function [row,col] = posAllImpl(obj)
            row=repmat((1:obj.shape(1))',obj.shape(2),1);
            col=reshape(repmat((1:obj.shape(2)),obj.shape(1),1),[],1);
        end
        function [nz] = nzPosImpl(obj,row,col)
            nz=sub2ind(obj.shape,row,col);
        end
    end
end