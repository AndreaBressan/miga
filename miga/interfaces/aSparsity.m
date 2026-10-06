classdef aSparsity <matlab.mixin.Heterogeneous
    % aSparsity is the interface that allows a sparse representation
    %  in which the non-zero data are stored in a vector
    %  and their positions are stored or otherwise known

    properties(GetAccess=public,SetAccess=protected)
        % 'nnz' is the number of non-zeros
        %
        % 'shape' is the size of the sparsity pattern
        %

        nnz     (1,1) double
        shape   (1,2) double
    end

    methods(Sealed)
        function [dir] = defaultDir(obj,data)
            %DEFAULTDIR Determine the default direction containing the basis data.
            %
            %   DIR = obj.defaultDir(DATA)
            %   determines the dimension of DATA along which the entries associated
            %   with the nonzero elements of OBJ are stored.
            %
            %   The method uses the number of nonzero entries OBJ.nnz and the
            %   dimensions of DATA to identify the appropriate direction. In the
            %   standard case, the last dimension of DATA is assumed to contain the
            %   OBJ.nnz entries.
            %
            %   If OBJ contains a single nonzero entry and the last dimension of
            %   DATA contains more than one element, the direction is shifted by one
            %   dimension. Conversely, if the last dimension contains a single
            %   element while OBJ.nnz is greater than one, the direction is shifted
            %   one dimension backwards.
            %
            %   Input arguments:
            %       OBJ  - aSparsity object containing the number of nonzero entries, OBJ.nnz.
            %
            %       DATA - array whose dimensions are inspected to determine the
            %              direction containing the OBJ.nnz entries.
            %
            %   Output arguments:
            %       DIR  - dimension of DATA expected to contain OBJ.nnz entries.


            dir=ndims(data);
            if obj.nnz==1 && size(data,dir)>1; dir=dir+1; end
            if size(data,dir)==1 && obj.nnz>1; dir=dir-1; end
            assert(size(data,dir)==obj.nnz, "Data has wrong size in the last direction")
        end

        function [out] = toFull(obj,data,dir)
            %TOFULL Expand data from sparse/compressed to full representation.
            %
            %   OUT = obj.toFull(DATA)
            %   expands DATA from the compressed representation associated with OBJ
            %   to its full representation.
            %
            %   OUT = obj.toFull(DATA,DIR)
            %   specifies the dimension DIR of DATA along which the compressed
            %   representation is stored.
            %
            %   If DATA has dimensions
            %
            %       (n1,...,n_(DIR-1),NNZ,n_(DIR+1),...,nM),
            %
            %   where NNZ = obj.nnz, the corresponding dimension in OUT is replaced
            %   by the dimensions of the sparsity structure represented by OBJ.
            %   The entries are placed in their corresponding positions, with zeros
            %   inserted where no active entry is present.
            %
            %   Input arguments:
            %       OBJ  - aSparsity object describing the sparsity pattern used for
            %              the compressed representation.
            %
            %       DATA - array containing the data in compressed form. Its size
            %              along dimension DIR must be equal to OBJ.nnz.
            %
            %       DIR  - dimension of DATA along which the compressed entries are
            %              stored. If omitted, the dimension is determined
            %              automatically by OBJ.defaultDir(DATA).
            %
            %   Output arguments:
            %       OUT  - data expanded to the full representation described by
            %              OBJ. The entries corresponding to the nonzero positions of
            %              OBJ contain the values from DATA, while the remaining
            %              entries are set to zero.

            arguments
                obj (1,1) aSparsity
                data      double
                dir (1,1) double = obj.defaultDir(data)
            end
            assert(size(data,dir)==obj.nnz, "Incorrect data size in the specified direction")
            out=toFullImpl(obj,data,dir);
        end

        function [out] = toSparse(obj,data,dir)
            %TOSPARSE Convert compressed data to MATLAB sparse matrix format.
            %
            %   OUT = obj.toSparse(DATA)
            %   converts DATA from the compressed representation associated with OBJ
            %   into a MATLAB sparse matrix.
            %
            %   OUT = obj.toSparse(DATA,DIR)
            %   specifies the dimension DIR of DATA containing the compressed entries.
            %
            %   If DATA has dimensions
            %
            %       (n1,...,n_(DIR-1),NNZ),
            %
            %   where NNZ = obj.nnz, the first DIR-1 dimensions are flattened and
            %   the compressed data are mapped to the corresponding positions of the
            %   sparsity pattern represented by OBJ.
            %
            %   The resulting sparse matrix has size
            %
            %       prod(n1,...,n_(DIR-1)) * shape(1)  by  shape(2),
            %
            %   where shape is the size of the sparsity pattern represented by OBJ.
            %
            %   The conversion is performed with the convention that DIR is the last
            %   non-singleton dimension of DATA. Dimensions after DIR must therefore
            %   be singleton dimensions.
            %
            %   Input arguments:
            %       OBJ  - aSparsity object describing the sparsity pattern used for
            %              the compressed representation.
            %
            %       DATA - array containing the data in compressed form. The size of
            %              DATA along dimension DIR must be equal to OBJ.nnz.
            %
            %       DIR  - dimension of DATA containing the compressed entries.
            %              If omitted, the dimension is determined automatically
            %              using OBJ.defaultDir(DATA).
            %
            %   Output arguments:
            %       OUT  - MATLAB sparse matrix containing the data expanded according
            %              to the sparsity pattern OBJ. The dimensions preceding DIR
            %              are flattened into the row dimension.

            arguments
                obj (1,1)  aSparsity
                data       double
                dir (1,1)  double = obj.defaultDir(data)
            end
            assert(all(size(data,dir+1:ndims(data))==1), "toSparse only works with the last dimension")
            [out] = toSparseImpl(obj,data,dir);
        end

        function [out] = toCompressed(obj,data,dir)
            %TOCOMPRESSED Convert data from full to compressed representation.
            %
            %   OUT = obj.toCompressed(DATA,DIR)
            %   extracts from DATA the entries corresponding to the nonzero positions
            %   of the sparsity pattern represented by OBJ and stores them in a
            %   compressed representation along dimension DIR.
            %
            %   If DATA has dimensions
            %
            %       (n1,...,nM),
            %
            %   the two dimensions DIR and DIR+1 must correspond to the dimensions
            %   of the sparsity pattern OBJ. The output has dimensions
            %
            %       (n1,...,n_(DIR-1),NNZ,n_(DIR+2),...,nM),
            %
            %   where NNZ = obj.nnz. The entries corresponding to the nonzero
            %   positions of the sparsity pattern are extracted and stored along
            %   dimension DIR.
            %
            %   Input arguments:
            %       OBJ  - aSparsity object describing the sparsity pattern used to
            %              identify the entries to be extracted.
            %
            %       DATA - full array containing the data to be compressed. The
            %              dimensions DIR and DIR+1 must have size equal to
            %              OBJ.shape.
            %
            %       DIR  - dimension along which the compressed representation is
            %              stored. The dimensions DIR and DIR+1 of DATA correspond
            %              to the two dimensions of the sparsity pattern.
            %
            %   Output arguments:
            %       OUT  - data in compressed representation. The two dimensions
            %              corresponding to the sparsity pattern are replaced by a
            %              single dimension of size OBJ.nnz.
            %
            %   The function checks that the size of DATA along dimensions DIR and
            %   DIR+1 is compatible with the sparsity pattern before performing the
            %   conversion.

            arguments
                obj (1,1)  aSparsity
                data       double
                dir (1,1)  double
            end
            assert(isequal(size(data,[dir,dir+1]),obj.shape),"toCompressed: shape mismatch between sparsity and data argument")
            [out] = toCompressedImpl(obj,data,dir);
        end

        function [out] = nnzCol(obj,col)
            %NNZCOL Number of nonzero entries in selected columns.
            %
            %   OUT = obj.nnzCol(COL)
            %   returns the number of nonzero entries in each column specified by COL.
            %   If COL contains multiple column indices, OUT is a column vector whose
            %   entries correspond to the selected columns.
            %
            %   OUT = obj.nnzCol()
            %   returns the number of nonzero entries in every column of the sparsity
            %   pattern represented by OBJ.
            %
            %   Input arguments:
            %       COL - Column indices for which the number of nonzero entries is
            %             requested. COL must contain positive integers. If omitted or
            %             empty, all columns are considered.
            %
            %   Output arguments:
            %       OUT - Column vector containing the number of nonzero entries in the
            %             requested columns. Its length is equal to numel(COL) when COL
            %             is specified, or to the number of columns of the sparsity
            %             pattern when COL is empty.


            arguments(Input)
                obj (1,1) aSparsity
                col {mustBeInteger,mustBePositive} = []
            end
            arguments(Output)
                out (:,1) double
            end
            if isempty(col); out=obj.nnzColAllImpl();
            else, out=obj.nnzColImpl(col); end
        end

        function [out] = nnzRow(obj,row)
            %NNZROW Number of nonzero entries in selected rows.
            %
            %   OUT = obj.nnzRow(ROW)
            %   returns the number of nonzero entries in each row specified by ROW.
            %   ROW can contain one or more positive integer row indices. OUT is a
            %   column vector with one entry for each specified row.
            %
            %   OUT = obj.nnzRow()
            %   returns the number of nonzero entries in every row of the sparsity
            %   pattern represented by OBJ.
            %
            %   Input arguments:
            %       ROW - Row indices for which the number of nonzero entries is
            %             requested. ROW must contain positive integers. If ROW is
            %             omitted or empty, all rows are considered.
            %
            %   Output arguments:
            %       OUT - Column vector containing the number of nonzero entries in the
            %             selected rows, or in all rows if ROW is empty.

            arguments(Input)
                obj (1,1) aSparsity
                row {mustBeInteger,mustBePositive} = []
            end
            arguments(Output)
                out (:,1) double
            end
            if isempty(row); out=obj.nnzRowAllImpl();
            else, out=obj.nnzRowImpl(row); end
        end

        function [row,col] = pos(obj,nz)
            %POS Row and column positions of selected nonzero entries.
            %
            %   [ROW,COL] = obj.pos(NZ)
            %   returns the row and column positions of the nonzero entries specified
            %   by NZ. If NZ contains multiple indices, ROW and COL are column vectors
            %   containing the corresponding row and column positions.
            %
            %   [ROW,COL] = obj.pos()
            %   returns the row and column positions of all nonzero entries of the
            %   sparsity pattern represented by OBJ.
            %
            %   Input arguments:
            %       NZ  - Indices of the nonzero entries whose positions are requested.
            %             NZ must contain positive integers. If omitted or empty, the
            %             positions of all nonzero entries are returned.
            %
            %   Output arguments:
            %       ROW - Column vector containing the row indices of the selected
            %             nonzero entries.
            %
            %       COL - Column vector containing the column indices of the selected
            %             nonzero entries.

            arguments(Input)
                obj (1,1) aSparsity
                nz {mustBeInteger,mustBePositive} = []
            end
            arguments(Output)
                row(:,1) double
                col(:,1) double
            end
            if isempty(nz); [row,col]=obj.posAllImpl();
            else, [row,col]=obj.posImpl(nz(:));end
        end

        function [nz] = nzPos(obj,row,col)
            %NZPOS Positions in the data vector corresponding to row and column indices.
            %
            %   NZ = obj.nzPos(ROW,COL)
            %   returns the positions in the data vector corresponding to the specified
            %   row and column indices of the sparsity pattern represented by OBJ.
            %   ROW and COL must have the same size, and each pair (ROW,COL) identifies
            %   one matrix entry.
            %
            %   NZ = obj.nzPos()
            %   returns the positions corresponding to the default empty ROW and COL
            %   inputs, as determined by the underlying implementation.
            %
            %   Input arguments:
            %       OBJ  - aSparsity object.
            %
            %       ROW - Row indices of the entries whose positions are requested.
            %             ROW must contain positive integers.
            %
            %       COL - Column indices of the entries whose positions are requested.
            %             COL must contain positive integers and must have the same
            %             size as ROW.
            %
            %   Output arguments:
            %       NZ  - Column vector containing the positions in the data vector
            %             associated with the specified row-column pairs.

            arguments(Input)
                obj (1,1) aSparsity
                row {mustBeInteger,mustBePositive} = []
                col {mustBeInteger,mustBePositive} = []
            end
            arguments(Output)
                nz (:,1) double
            end
            assert(isequal(size(row),size(col)),"nzPos, row and col arguments must have the same shape")
            nz=obj.nzPosImpl(row,col);
        end
    end

    methods
        function out=toFullImpl(obj,data,dir)
            %TOFULLIMPL Expand compressed data to its full representation.
            %
            %   OUT = obj.toFullImpl(DATA,DIR)
            %   expands the compressed dimension of DATA according to the sparsity
            %   pattern represented by OBJ. The entries of DATA are placed in the
            %   corresponding row and column positions of the full representation,
            %   while all remaining entries are set to zero.
            %
            %   If DATA has size
            %
            %       [n1,...,n_(DIR-1),NNZ,n_(DIR+1),...,nM],
            %
            %   where NNZ = obj.nnz, then OUT has size
            %
            %       [n1,...,n_(DIR-1),obj.shape(1),obj.shape(2),
            %        n_(DIR+1),...,nM].
            %
            %   The NNZ entries along dimension DIR are mapped to the nonzero positions
            %   of the sparsity pattern OBJ.
            %
            %   Input arguments:
            %       OBJ  - aSparsity object.
            %
            %       DATA - Array containing the data in compressed representation, with
            %              obj.nnz entries along dimension DIR.
            %
            %       DIR  - Dimension along which the compressed data is stored.
            %
            %   Output arguments:
            %       OUT  - Full array obtained by expanding the compressed dimension
            %              according to the sparsity pattern OBJ. Entries not belonging
            %              to the sparsity pattern are set to zero.

            [row,col]=obj.posAllImpl();
            sz=size(data,1:max(ndims(data),dir));
            out=zeros([prod(sz(1:(dir-1))),obj.shape,prod(sz((dir+1):ndims(data)))]);
            dataR=reshape(data,[prod(sz(1:(dir-1))),obj.nnz,prod(sz((dir+1):ndims(data)))]);
            for i=1:numel(row)
                out(:,row(i),col(i),:)=dataR(:,i,:);
            end
            out=reshape(out,[sz(1:(dir-1)),obj.shape,sz((dir+1):ndims(data))]);
        end

        function out=toSparseImpl(obj,data,dir)
            %TOSPARSEIMPL Convert compressed data to a MATLAB sparse matrix.
            %
            %   OUT = obj.toSparseImpl(DATA,DIR)
            %   converts DATA from its compressed representation into a MATLAB sparse
            %   matrix using the sparsity pattern represented by OBJ.
            %
            %   If DATA has size
            %
            %       [n1,...,n_(DIR-1),NNZ,n_(DIR+1),...,nM],
            %
            %   where NNZ = obj.nnz, the dimensions preceding DIR are flattened into a
            %   single dimension. The resulting sparse matrix has size
            %
            %       [prod(n1,...,n_(DIR-1))*obj.shape(1), obj.shape(2)].
            %
            %   The entries of DATA are placed in the row and column positions
            %   corresponding to the nonzero entries of the sparsity pattern OBJ.
            %   When multiple entries precede dimension DIR, the corresponding row
            %   indices are replicated with the appropriate offset.
            %
            %   Input arguments:
            %       OBJ  - aSparsity object.
            %
            %       DATA - Array containing the data in compressed representation, with
            %              obj.nnz entries along dimension DIR.
            %
            %       DIR  - Dimension along which the compressed data is stored.
            %
            %   Output arguments:
            %       OUT  - MATLAB sparse matrix containing the data in the positions
            %              specified by the sparsity pattern OBJ.

            sz=size(data,1:max(ndims(data),dir));
            [row,col]=obj.posAllImpl();
            rep=prod(sz(1:(dir-1)));
            if rep>1
                row=(1:rep)'+rep*(row(:)'-1);
                col=repelem(col(:),rep);
            end
            out=sparse(row(:),col,data(:),obj.shape(1)*rep,obj.shape(2));
        end

        function out=toCompressedImpl(obj,data,dir)
            %TOCOMPRESSEDIMPL Convert full data to compressed representation.
            %
            %   OUT = obj.toCompressedImpl(DATA,DIR)
            %   extracts from DATA the entries corresponding to the nonzero positions
            %   of the sparsity pattern represented by OBJ and stores them in a
            %   compressed representation.
            %
            %   If DATA has size
            %
            %       [n1,...,n_(DIR-1),m1,m2,n_(DIR+2),...,nM],
            %
            %   where [m1,m2] = obj.shape, then OUT has size
            %
            %       [n1,...,n_(DIR-1),obj.nnz,n_(DIR+2),...,nM].
            %
            %   The two dimensions DIR and DIR+1 of DATA are therefore replaced by a
            %   single dimension of size obj.nnz. The entries along this dimension are
            %   ordered according to the nonzero positions of the sparsity pattern OBJ.
            %
            %   Input arguments:
            %       OBJ  - aSparsity object.
            %
            %       DATA - Array containing the full data representation. Dimensions
            %              DIR and DIR+1 must have sizes given by obj.shape.
            %
            %       DIR  - First of the two dimensions corresponding to the sparsity
            %              pattern.
            %
            %   Output arguments:
            %       OUT  - Array containing only the entries corresponding to the
            %              nonzero positions of OBJ, with the two sparsity dimensions
            %              compressed into a single dimension of size obj.nnz.

            dataSz=size(data,1:max(dir+1,ndims(data)));
            outSz=[dataSz(1:dir-1),obj.nnz,dataSz(dir+2:end),1];
            data=reshape(data,[prod(dataSz(1:dir-1)),dataSz(dir),dataSz(dir+1),prod(dataSz(dir+2:end))]);
            out=zeros([size(data,1),obj.nnz,size(data,4)]);
            [row,col]=obj.posAllImpl();
            % TODO replace loop with smarter thing
            for nz=1:obj.nnz
                out(:,nz,:)=data(:,row(nz),col(nz),:);
            end
            out=reshape(out,outSz);
        end

        function [out] = nnzColAllImpl(obj)
            out=obj.nnzColImpl(1:obj.shape(2));
        end

        function [out] = nnzRowAllImpl(obj)
            out=obj.nnzRowImpl(1:obj.shape(1));
        end

        function [row,col] = posAllImpl(obj)
            [row,col]=obj.posImpl(1:obj.nnz);
        end
    end

    methods(Abstract)
        [out] = nnzColImpl(obj,col)
        [out] = nnzRowImpl(obj,row)
        [out] = posImpl(obj,nnz)
        [nz]  = nzPosImpl(obj,row,col)
    end
end


