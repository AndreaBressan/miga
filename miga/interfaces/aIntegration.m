classdef aIntegration < matlab.mixin.Heterogeneous
    % aIntegration is the interface for quadratures and integration
    % methods.

    properties(GetAccess=public,SetAccess=protected)
        % 'domain' is a n x 2 matrix describing a n-dimensional
        %          axis aligned cuboid in R^n. The i-th row represents the
        %          i-th axis, while the columns are respectively the vertex
        %          having minimal and maximal coordinate.
        %
        % 'nodes' describe the qudrature nodes
        %
        % 'weights' contain the quadrature weights

        domain (:,2) double
        nodes
        weights
    end

    methods
        function [val,sparsity]=integrate(obj,weight,test,request,trial)
            %INTEGRATE Integrate weighted basis functions and their derivatives.
            %
            %   VAL = obj.integrate(WEIGHT)
            %   integrates the functions in WEIGHT over the domain and returns the
            %   resulting values.
            %
            %   VAL = obj.integrate(WEIGHT,TEST)
            %   integrates the product of WEIGHT and TEST over the domain and returns
            %   the resulting values.
            %
            %   VAL = obj.integrate(WEIGHT,TEST,REQUEST)
            %   integrates WEIGHT multiplied by the partial derivatives of TEST
            %   specified by REQUEST. The contributions associated with the different
            %   requested derivatives are summed, and the result is returned as a
            %   vector.
            %
            %   VAL = obj.integrate(WEIGHT,TEST,REQUEST,TRIAL)
            %   integrates WEIGHT multiplied by the requested partial derivatives of
            %   TEST and TRIAL. The result is a matrix whose entries correspond to the
            %   interactions between the basis functions of TEST and TRIAL.
            %
            %   If REQUEST is empty, the values of the basis functions themselves are
            %   used. Thus, in the matrix case, the default operation corresponds to
            %   the integral of WEIGHT times TEST times TRIAL.
            %
            %   [VAL,SPARSITY] = obj.integrate(...)
            %   also returns the sparsity pattern associated with the computed values.
            %   In the matrix case, SPARSITY describes the positions of the nonzero
            %   entries of the Galerkin matrix.
            %
            %   Input arguments:
            %       WEIGHT  - Set of weighting functions, represented by an object
            %                 deriving from aBasis.
            %
            %       TEST    - Scalar-valued basis used as test functions. TEST must
            %                 contain at most one basis object from aBasis.
            %
            %       REQUEST - Array of size [domainDim,numTerms,s], where domainDim is
            %                 the dimension of the computational domain. Each column
            %                 REQUEST(:,i,k) specifies the multi-index of the partial
            %                 derivative used in term i. If TRIAL is empty, s must be
            %                 equal to 1; otherwise, s must be equal to 2. In the
            %                 matrix case, REQUEST(:,i,1) specifies the derivative of
            %                 the test function and REQUEST(:,i,2) the derivative of
            %                 the trial function. The contributions associated with
            %                 different terms i are summed.
            %
            %       TRIAL   - Scalar-valued basis used as trial functions. TRIAL must
            %                 contain at most one basis object from aBasis.
            %
            %   Output arguments:
            %       VAL      - Result of the integration. Its structure depends on the
            %                  supplied arguments. When TEST and TRIAL are both
            %                  provided, VAL contains the entries of the resulting
            %                  Galerkin matrix. If only WEIGHT and TEST are provided,
            %                  it contains the corresponding vector of integrals.
            %
            %       SPARSITY - aSparsity object describing the positions of the
            %                  nonzero entries represented by VAL. In the matrix case,
            %                  it corresponds to the sparsity pattern between TEST and
            %                  TRIAL.
            %
            %   Notes:
            %       WEIGHT must be a basis object. TEST and TRIAL, when provided, must
            %       be scalar-valued bases.
            %
            %       The number of derivative terms specified by REQUEST must match the
            %       dimension of WEIGHT. In the matrix case, WEIGHT must be scalar-
            %       valued, i.e., its rangeSize must be [] or [1].
            %
            %       When only WEIGHT is provided, the result is a scalar. When TEST is
            %       provided without TRIAL, the result is a vector. When both TEST and
            %       TRIAL are provided, the result represents a matrix and, for a
            %       scalar-valued WEIGHT, is returned as a MATLAB sparse matrix when
            %       only VAL is requested.


            arguments
                obj
                weight   (1,1)   aBasis
                test     (:,:)   aBasis  =aBasis.empty
                request  (:,:,:) double  =double.empty
                trial    (:,:)   aBasis  =aBasis.empty
            end
            assert(all([numel(test), numel(trial)]<=1),'test and trial cannot be arrays')
            if ~isempty(request)
                assert(~isempty(test), "Cannot request partial derivatives of test when test is empty")
                assert(size(request,2)==weight.dimension, "The number of request terms should equal the dimension of the weight space")
                if isempty(trial)
                    assert( size(request,3)==1, "Request should contain a single page when trial is empty");
                else
                    assert( size(request,3)==2, "Request should contain a page for test and one for trial");
                end
            end
            assert(isempty(test)  || prod( test.rangeSize)==1,  "test must be a scalar valued basis");
            assert(isempty(trial) || prod(trial.rangeSize)==1, "trial must be a scalar valued basis");

            dim=size(obj.domain,1);
            wShape=weight.rangeSize;

            if isempty(test)
                iop=obj.getIntegrateOperator();
                if nargout==2; sparsity=fSparsity([1,1]); end
                oShape=wShape;
            elseif isempty(trial)
                if isempty(request); request=zeros(dim,1,1); end
                iop=obj.getTestOperator(test,request);
                oShape=[wShape, test.dimension,1];
                if nargout==2; sparsity=fSparsity([test.dimension,1]); end
            else
                if isempty(request); request=zeros(dim,1,2); end
                sparsity=test.sparsity(trial,obj.domain);
                [iop]=obj.getMatrixOperator(test,trial,request,sparsity);
                oShape=[wShape,sparsity.nnz];
            end

            wVals=weight.evaluate(obj.nodes);
            wVals=reshape(wVals,[prod(wShape),npts(obj.nodes),weight.dimension]);
            val=aIntegration.applyIOP(wVals,iop);
            val=reshape(val,[oShape,numel(val)/prod(oShape),1]);

            if nargout==1 && ~isempty(trial) && prod(wShape)==1
                val=sparsity.toSparse(val(:),1);
            end
        end

        function [fun]=getApplyFun(obj,weight,test,request,trial,options)
            %GETAPPLYFUN Construct a function handle for applying an assembled operator.
            %
            %   FUN = obj.getApplyFun(WEIGHT,TEST,REQUEST,TRIAL)
            %   constructs a function handle that maps a vector of coefficients in the
            %   TRIAL space to the corresponding vector in the TEST space. The
            %   resulting function performs the same operation as assembling the
            %   matrix associated with the bilinear form and multiplying it by the
            %   input vector, without explicitly assembling the matrix.
            %
            %   The operator is defined by the weighting functions WEIGHT and by the
            %   partial derivatives of TEST and TRIAL specified by REQUEST.
            %
            %   FUN = obj.getApplyFun(...,OPTIONS)
            %   allows additional options to be specified. In particular,
            %
            %       OPTIONS.cacheWeight = true
            %
            %   precomputes the values of WEIGHT at the integration nodes and reuses
            %   them every time FUN is evaluated. If false, the weight is evaluated
            %   each time FUN is called.
            %
            %   Input arguments:
            %       WEIGHT  - Scalar-valued basis defining the weighting function.
            %
            %       TEST    - Scalar-valued test basis.
            %
            %       REQUEST - Array specifying the partial derivatives used in the
            %                 operator. The derivatives of the test and trial
            %                 functions are specified by REQUEST(:,:,1) and
            %                 REQUEST(:,:,2), respectively.
            %
            %       TRIAL   - Scalar-valued trial basis. The input of the function
            %                 handle FUN is a vector of coefficients in this basis.
            %
            %       OPTIONS - Structure containing additional options.
            %
            %           cacheWeight - Logical flag indicating whether the values of
            %                         WEIGHT at the integration nodes should be
            %                         precomputed and cached. Default: true.
            %
            %   Output arguments:
            %       FUN - Function handle that takes a vector of trial coefficients
            %             and returns the corresponding vector in the test space.
            %             The action is equivalent to applying the matrix obtained by
            %             assembling the corresponding bilinear form.
            %
            %   Notes:
            %       The method is implemented only for scalar problems. WEIGHT, TEST,
            %       and TRIAL must therefore all have scalar range.
            %
            %       Caching the weight values can reduce the cost of repeated
            %       applications of FUN when the same weighting function is used
            %       multiple times.


            arguments
                obj
                weight   (1,1)   aBasis
                test     (:,:)   aBasis  =aBasis.empty
                request  (:,:,:) double  =double.empty
                trial    (:,:)   aBasis  =aBasis.empty
                options.cacheWeight (1,1) logical =true
            end
            assert(all([prod(weight.rangeSize),prod(test.rangeSize),prod(trial.rangeSize)]==1),"Apply is only implemented for scalar problems");
            if options.cacheWeight
                [w]=evaluateWeights(obj,weight);
                nW=@(tr_co)reshape(trial.evaluateCombination(tr_co,obj.nodes,parDers(request(:,:,2)))',1,[],size(request,2)).*w;
            else
                nW=@(tr_co)reshape(trial.evaluateCombination(tr_co,obj.nodes,parDers(request(:,:,2)))',1,[],size(request,2)).*evaluateWeights(obj,weight);
            end
            [iop]=obj.getApplyOperator(test,trial,request);
            fun=@(tr_co) aIntegration.applyIOP(nW(tr_co),iop);
        end
    end

    methods

        function [iOp]=getIntegrateOperator(obj)
            %GETINTEGRATEOPERATOR Return the integration weights.
            %
            %   IOP = obj.getIntegrateOperator()
            %   returns the integration weights stored in OBJ. The returned object is
            %   used as the integration operator when integrating a weighting function
            %   without test or trial functions.
            %
            %   Output arguments:
            %       IOP - Integration weights associated with the quadrature rule.

            iOp=obj.weights;
        end

        function [iOp]=getTestOperator  (obj,test,request)
            %GETTESTOPERATOR Construct the integration operator for test functions.
            %
            %   IOP = obj.getTestOperator(TEST,REQUEST)
            %   constructs the integration operator associated with the test basis
            %   TEST and the derivative requests specified by REQUEST.
            %
            %   The method evaluates the requested derivatives of TEST at the
            %   integration nodes, converts the resulting basis evaluations to sparse
            %   form according to their active-point structure, and multiplies them by
            %   the quadrature weights.
            %
            %   Input arguments:
            %       TEST    - Basis of test functions to be evaluated at the
            %                 integration nodes.
            %
            %       REQUEST - Derivative requests passed to TEST.evaluateI. Each
            %                 requested derivative defines one component of the
            %                 resulting integration operator.
            %
            %   Output arguments:
            %       IOP - Cell array containing the integration operators associated
            %             with the requested derivatives. Each element combines the
            %             values of the corresponding derivative of TEST with the
            %             quadrature weights.
            %
            %   Notes:
            %       The quadrature weights are obtained from OBJ.weights. If the
            %       weights are stored as a cell array, the entries are combined
            %       componentwise to form the multidimensional quadrature weights.
            %
            %       The dimension of the evaluated test functions is taken into account
            %       when arranging the sparse evaluation data. If TEST has range size
            %       R, the resulting operators use blocks of size
            %       prod(R)*nPts, where nPts is the number of active integration
            %       points.

            dimRange=prod(test.rangeSize);
            [v,a]=test.evaluateI(obj.nodes,request);
            nPts=a.shape(1);
            nFun=a.shape(2);
            if ~iscell(obj.weights)
                wq=obj.weights;
            else
                wq=fold(@(a,b) a.*reshape(b,[ones(1,ndims(a)),numel(b)]),obj.weights,1);
            end
            wq=repelem(wq(:),dimRange,1);

            iOp=cellfun(@(w)wq.*reshape(a.toSparse(w),dimRange*nPts,nFun),v,"UniformOutput",false);
        end


        function [iop]=getApplyOperator(obj,test,trial,request)
            %GETAPPLYOPERATOR Construct the operator used for matrix-free application.
            %
            %   IOP = obj.getApplyOperator(TEST,TRIAL,REQUEST)
            %   constructs the integration operator used by getApplyFun to apply the
            %   operator associated with the test and trial bases without explicitly
            %   assembling the corresponding matrix.
            %
            %   The operator is constructed from the derivatives of the TEST basis
            %   specified by REQUEST(:,:,1). The evaluation of the TRIAL basis and its
            %   derivatives is incorporated separately when the resulting function
            %   handle is applied.
            %
            %   Input arguments:
            %       TEST    - Basis of test functions.
            %
            %       TRIAL   - Basis of trial functions. The argument is used to define
            %                 the trial space of the operator.
            %
            %       REQUEST - Array specifying the partial derivatives involved in the
            %                 operator. REQUEST(:,:,1) contains the derivative
            %                 requests associated with the test functions.
            %
            %   Output arguments:
            %       IOP - Integration operator associated with the requested
            %             derivatives of the test basis.

            iop=obj.getTestOperator(test,request(:,:,1));
        end


        function [iOp]=getMatrixOperator(obj,test,trial,request,sparsity)
            %GETMATRIXOPERATOR Construct sparse integration operators for a test-trial pair.
            %
            %   IOP = obj.getMatrixOperator(TEST,TRIAL,REQUEST,SPARSITY)
            %   constructs the sparse integration operators associated with the
            %   interactions between the TEST and TRIAL bases. The requested
            %   derivatives of the two bases are evaluated at the integration nodes,
            %   combined with the quadrature weights, and assembled according to the
            %   sparsity pattern SPARSITY.
            %
            %   The resulting operators contain one row for each integration point and
            %   one column for each nonzero entry of the test-trial sparsity pattern.
            %   Different requested derivative terms are returned as separate elements
            %   of the output cell array.
            %
            %   Input arguments:
            %       TEST     - Scalar-valued basis of test functions.
            %
            %       TRIAL    - Scalar-valued basis of trial functions.
            %
            %       REQUEST  - Three-dimensional array specifying the derivatives used
            %                  in the bilinear form. REQUEST(:,:,1) contains the
            %                  derivative requests for TEST, while REQUEST(:,:,2)
            %                  contains those for TRIAL.
            %
            %       SPARSITY - aSparsity object describing the nonzero interactions
            %                  between TEST and TRIAL.
            %
            %   Output arguments:
            %       IOP      - Cell array of sparse integration operators. Each element
            %                  corresponds to one derivative term in REQUEST and has
            %                  size [nPts,SPARSITY.nnz], where nPts is the number of
            %                  integration points.
            %
            %   Notes:
            %       The method is implemented only for scalar-valued TEST and TRIAL
            %       bases.
            %
            %       The columns of the operators are ordered according to the nonzero
            %       positions described by SPARSITY. The quadrature weights stored in
            %       OBJ are incorporated into the operators.

            assert(prod(test.rangeSize)==1 && prod(trial.rangeSize)==1,...
                "Only works with scalar bases")

            [vr,ar]= test.evaluateI(obj.nodes,request(:,:,1));
            [vc,ac]=trial.evaluateI(obj.nodes,request(:,:,2));

            nPts=ar.shape(1);

            [~,fn_r]=ar.pos();
            [~,fn_c]=ac.pos();
            pt_nnz_r=ar.nnzRow;
            pt_nnz_c=ac.nnzRow;
            pt_nnz=pt_nnz_r.*pt_nnz_c;
            pt_block_start_r=cumsum([1;pt_nnz_r]);
            pt_block_start=cumsum([1;pt_nnz]);

            [~,~,row_tmp]=find(ar.toSparse(1:ar.nnz)');
            row=zeros(pt_block_start(end)-1,1);
            for pt=1:nPts
                rs=row_tmp(pt_block_start_r(pt):pt_block_start_r(pt+1)-1);
                row(pt_block_start(pt):pt_block_start(pt+1)-1)=reshape(repmat((rs),1,pt_nnz_c(pt)),[],1);
            end
            [~,~,col_tmp]=find(ac.toSparse(1:ac.nnz)');
            col=reshape(repelem(col_tmp,repelem(pt_nnz_r,pt_nnz_c)),size(row));
            pts=repelem(1:nPts,pt_nnz);
            nnz=sparsity.nzPos(fn_r(row),fn_c(col));

            if ~iscell(obj.weights)
                qw=obj.weights(pts);
            else
                pt_c=cell(numel(obj.weights),1);
                [pt_c{:}]=ind2sub([cellfun(@numel,obj.weights),1,1], pts);
                tmp=cellfun(@(a,b) a(b), obj.weights,pt_c,"UniformOutput",false);
                qw=fold(@(a,b) a.*b, tmp, ones(nPts,1));
            end

            iOp=cellfun(@(r,c)sparse(pts,nnz,reshape(r(row),[],1).*reshape(c(col),[],1).*qw(:),nPts,sparsity.nnz),vr,vc,"UniformOutput",false);
        end
    end

    methods(Static)
        function val = applyIOP(m,iop)
            %APPLYIOP Apply an integration operator to input data.
            %
            %   VAL = applyIOP(M,IOP)
            %   applies the integration operator IOP to the data M. IOP can be
            %   represented by a dense or sparse matrix, by a cell array of matrices,
            %   or by an operator that can be applied directly to M.
            %
            %   If IOP is a dense matrix, the operation is performed using page-wise
            %   matrix multiplication. If IOP is sparse, the multiplication is
            %   performed separately on each page of M.
            %
            %   If IOP is a cell array, the matrices in the cell array are first
            %   converted into an operator and the resulting operator is applied to M.
            %
            %   If IOP is neither a numeric matrix nor a cell array, it is assumed to
            %   be an operator that accepts M as input, such as a function handle or
            %   a specialized operator for tensor-product structures.
            %
            %   Input arguments:
            %       M   - Input data to which the integration operator is applied.
            %
            %       IOP - Integration operator. It can be:
            %             * a dense or sparse matrix;
            %             * a cell array of matrices;
            %             * a function handle or another operator supporting direct
            %               application to M.
            %
            %   Output arguments:
            %       VAL - Result of applying IOP to M. Its dimensions depend on the
            %             representation of IOP and the dimensions of M.

            if isa(iop,'double')
                if ~issparse(iop)
                    val=pagemtimes(m,iop);
                else
                    val=zeros(size(m,1), size(iop,2), prod(size(m,3:ndims(m))));
                    for p=1:size(val,3)
                        val(:,:,p)=m(:,:,p)*iop;
                    end
                    val=reshape(val,[size(m,1), size(iop,2),size(m,3:ndims(m))]);
                end
            elseif iscell(iop)
                iOp=aIntegration.matricesToOp(iop);
                val=iOp(m);
            else
                val=iop(m); % special operators class, e.g. for tensor product
            end
        end


        function op=matricesToOp(matrices)
            %MATRICESTOOP Construct an operator from a cell array of matrices.
            %
            %   OP = matricesToOp(MATRICES)
            %   constructs a function handle that applies the sum of the operators
            %   represented by the matrices in MATRICES.
            %
            %   Each matrix MATRICES{I} acts on the I-th page of the input M:
            %
            %       M(:,:,I) * MATRICES{I}.
            %
            %   The resulting operator therefore computes
            %
            %       OP(M) = sum_I M(:,:,I) * MATRICES{I}.
            %
            %   The individual matrix applications are combined recursively in a
            %   balanced binary tree. This provides a compact representation of the
            %   resulting operator and avoids explicitly forming the sum.
            %
            %   Input arguments:
            %       MATRICES - Cell array containing the matrices defining the
            %                  individual page-wise operators.
            %
            %   Output arguments:
            %       OP       - Function handle that takes an array M as input and
            %                  returns the sum of the matrix products associated with
            %                  the pages of M.


            n=numel(matrices);
            tmp=arrayfun(@(i)@(m) m(:,:,i)*matrices{i},1:n,"UniformOutput",false);
            while n>1
                mid=floor(n/2);
                for i=1:mid
                    tmp{i}=@(m) tmp{2*i-1}(m)+tmp{2*i}(m);
                end
                if mod(n,2)
                    tmp{mid+1}=tmp{end};
                    n=mid+1;
                else
                    n=mid;
                end
            end
            op=tmp{1};
        end
    end
end



function [w]=evaluateWeights(Q,weight)
%EVALUATEWEIGHTS Evaluate weighting functions at quadrature nodes.
%
%   W = evaluateWeights(Q,WEIGHT)
%   evaluates the weighting functions WEIGHT at the integration nodes
%   stored in Q and reshapes the result according to the range size of
%   WEIGHT and the number of integration nodes.
%
%   Input arguments:
%       Q      - Integration object containing the quadrature nodes in
%                Q.nodes.
%
%       WEIGHT - Basis defining the weighting functions to be evaluated.
%
%   Output arguments:
%       W      - Values of WEIGHT at the integration nodes. The result has
%                size
%
%                    [prod(WEIGHT.rangeSize),nPts(Q.nodes),WEIGHT.dimension].
%
%                The first dimension contains the components of the range,
%                the second indexes the integration nodes, and the third
%                indexes the basis functions.


wShape=weight.rangeSize;
w=weight.evaluate(Q.nodes);
w=reshape(w,[prod(wShape),npts(Q.nodes),weight.dimension]);
end