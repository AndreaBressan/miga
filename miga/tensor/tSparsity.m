
classdef tSparsity <aSparsity
    % tSparsity, represents the sparsity that one gets using tensor product
    % integration, i.e. sum factorization.
    properties
        uSparsity
        uShape
        uNnz
    end

    methods
        function obj=tSparsity(components)
            %TSPARSITY Construct a tensor-product sparsity pattern.
            %
            %   OBJ = tSparsity(COMPONENTS)
            %   constructs a tensor-product sparsity pattern from the component
            %   sparsity patterns specified in COMPONENTS.
            %
            %   The resulting sparsity pattern represents the tensor product of the
            %   sparsity patterns of the individual components.
            %
            %   Input arguments:
            %       COMPONENTS - Column vector of aSparsity objects defining the
            %                    component sparsity patterns. Each component
            %                    corresponds to one direction of the tensor-product
            %                    structure.
            %
            %   Output arguments:
            %       OBJ - Tensor-product sparsity object. Its main properties are:
            %
            %             uSparsity - Component sparsity patterns.
            %             uShape    - Shapes of the component sparsity patterns.
            %             shape     - Overall shape of the tensor-product sparsity
            %                         pattern.
            %             uNnz      - Number of nonzero entries of each component.
            %             nnz       - Total number of nonzero entries of the
            %                         tensor-product sparsity pattern.
            %
            %   Notes:
            %       If the component sparsity patterns have shapes
            %
            %           [m_i,n_i],  i = 1,...,D,
            %
            %       then the tensor-product sparsity pattern has shape
            %
            %           [prod_i m_i, prod_i n_i].
            %
            %       Similarly, if the i-th component if the i-th component contains nnz_i nonzero entries,
            % the tensor-product sparsity pattern contains prod_i nnz_i nonzero entries.


            arguments
                components (:,1) aSparsity
            end
            obj.uSparsity=components;
            obj.uShape=cat(1,components.shape);
            obj.shape=prod(obj.uShape,1);
            obj.uNnz=cat(1,components.nnz);
            obj.nnz=prod(obj.uNnz);
        end

        function [out] = nnzColImpl(obj,col)
            dim=numel(obj.uSparsity);
            ind=cell(dim,1);
            [ind{:}]=ind2sub(obj.uShape(:,2),col(:));
            out=1;
            for d=1:dim
                out=out.*obj.uSparsity(d).nnzColImpl(ind{d});
            end
        end

        function [out] = nnzColAllImpl(obj)
            dim=numel(obj.uSparsity);
            out=1;
            for d=1:dim
                univariate=obj.uSparsity(d).nnzColAllImpl();
                out=out(:).*univariate(:)';
            end
            out=out(:);
        end

        function [out] = nnzRowImpl(obj,row)
            dim=numel(obj.uSparsity);
            ind=cell(dim,1);
            [ind{:}]=ind2sub(obj.uShape(:,1),row(:));

            out=1;
            for d=1:dim
                out=out.*obj.uSparsity(d).nnzRowImpl(ind{d});
            end
        end

        function [out] = nnzRowAllImpl(obj)
            dim=numel(obj.uSparsity);
            out=1;
            for d=1:dim
                univariate=obj.uSparsity(d).nnzRowAllImpl();
                out=out(:).*univariate(:)';
            end
            out=out(:);
        end

        function [row,col] = posImpl(obj,nnz)
            dim=numel(obj.uSparsity);
            shiftR=cumprod([1,obj.uShape(:,1)']);
            shiftC=cumprod([1,obj.uShape(:,2)']);
            nnzMI=cell(dim,1);
            [nnzMI{:}]=ind2sub(obj.uNnz,nnz);
            row=1;
            col=1;
            for d=1:dim
                [rowComp, colComp]=obj.uSparsity(d).posImpl(nnzMI{d});
                row=row+(rowComp -1)*shiftR(d);
                col=col+(colComp -1)*shiftC(d);
            end
        end

        function [nz]=nzPosImpl(obj,row,col)
            dim=size(obj.uShape,1);
            uRow=cell(1,dim); [uRow{:}]=ind2sub(obj.uShape(:,1)',row);
            uCol=cell(1,dim); [uCol{:}]=ind2sub(obj.uShape(:,2)',col);
            uNz =arrayfun(@(d)obj.uSparsity(d).nzPosImpl(uRow{d},uCol{d}),1:dim,'UniformOutput',false );
            nz=reshape(sub2ind(obj.uNnz,uNz{:}),size(row));
        end

        function [row,col] = posAllImpl(obj)
            dim=numel(obj.uSparsity);
            shiftR=cumprod([1,obj.uShape(:,1)']);
            shiftC=cumprod([1,obj.uShape(:,2)']);
            row=1;
            col=1;
            for d=1:dim
                [rowComp, colComp]=obj.uSparsity(d).posAllImpl();
                row=row(:)+(rowComp(:)' -1)*shiftR(d);
                col=col(:)+(colComp(:)' -1)*shiftC(d);
            end
            row=row(:);
            col=col(:);
        end

        function out=toSparseImpl(obj,data,dir)
            if 3==exist('tSparsityToSparse','file') && numel(data)==obj.nnz % use mex file
                % TODO fix mex file working only for column major component sparsities
                dim=numel(obj.uSparsity);
                sparsityData=cell(3,dim);
                for i=1:dim
                    [sparsityData{1,i},sparsityData{2,i}]=obj.uSparsity(i).pos();
                    sparsityData{3,i}=obj.uSparsity(i).nnzCol();
                end
                out=tSparsityToSparse(sparsityData,obj.uShape, data, maxNumCompThreads);
            else % use generic implementation
                out=obj.toSparseImpl@aSparsity(data,dir);
            end
        end
    end
end