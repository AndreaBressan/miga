classdef kronMat
    properties (SetAccess = private)
        terms (1,:) cell {mustBeAllMatrices} = {}
        trans (1,1) logical = false
    end

    methods

        function obj = kronMat(varargin)
            if nargin == 1 && isa(varargin{1},'cell')
                obj.terms = varargin{1};
            else
                obj.terms = varargin;
            end
        end

        function obj = ctranspose(obj)
            obj.trans = ~obj.trans;
            obj.terms = cellfun(@conj,obj.terms,'UniformOutput',false);
        end

        function obj = transpose(obj)
            obj.trans = ~obj.trans;
        end

        function [varargout] = shape(obj)
            rows = cellfun(@(x) size(x,1), obj.terms');
            cols = cellfun(@(x) size(x,2), obj.terms');
            if obj.trans
                [cols,rows] = deal(rows,cols);
            end
            if nargout <= 1
                varargout{1} = [rows,cols];
            elseif nargout == 2
                varargout{1} = rows;
                varargout{2} = cols;
            else
                error('kronMat:OutputMismatch', ...
                    'Too many output arguments. "shape" returns either [rows, cols] or separate row and column vectors.');
            end
        end

        function varargout = size(obj,dim)
            arguments
                obj
                dim (1,:) {mustBeMember(dim, [1, 2])} = [1,2]
            end
            [r, c] = obj.shape();
            fullSize = [prod(r),prod(c)];
            if nargout <= 1
                varargout{1} = fullSize(dim);
            else
                assert(nargout == numel(dim),'kronMat:OutputMismatch',...
                    'Number of output arguments must equal the number of input dimension arguments.');
                varargout = arrayfun(@(x) fullSize(x),dim,'UniformOutput',false);
            end
        end

        function v = mtimes(obj,B)
            if isa(obj,'kronMat') && isa(B,'double')
                v = apply(obj,B);
            elseif isa(obj,'kronMat') && isa(B,'kronMat')
                assert(numel(obj.terms) == numel(B.terms), ...
                    'kronMat:TermMismatch','Number of terms must match.');
                n = numel(obj.terms);
                newTerms = cell(1, n);
                for i = 1:n
                    newTerms{i} = obj.terms{i}*B.terms{i};
                end
                v = kronMat(newTerms);
            else
                error(['Product between kronMat and ',class(B),' is not implemented.']);
            end
        end

        function v = mldivide(obj,B)
            if isa(obj,'kronMat') && isa(B,'double')
                v = applyInv(obj,B);
            elseif isa(obj,'kronMat') && isa(B,'kronMat')
                assert(numel(obj.terms) == numel(B.terms), ...
                    'kronMat:FactorMismatch','Number of factors must match.');
                n = numel(obj.terms);
                newTerms = cell(1, n);
                for i = 1:n
                    newTerms{i} = obj.terms{i}\B.terms{i};
                end
                v = kronMat(newTerms);
            else
                error(['Product between ',class(B),' and kronMat is not implemented.']);
            end
        end

        function obj = kronPower(obj,k)
            arguments
                obj kronMat
                k (1,1) {mustBeNumeric, mustBeInteger, mustBePositive}
            end
            assert(numel(obj.terms) == 1,'kronPower only works on a single base term.');
            obj.terms = repmat(obj.terms,1,k);
        end

        function M = toMat(obj)
            M = 1;
            for i = 1:numel(obj.terms)
                if obj.trans
                    term = obj.terms{i}.';
                else
                    term = obj.terms{i};
                end
                M = kron(M,term);
            end
        end

        function v = apply(obj,x)
            assert(obj.size(2) == size(x,1), 'kronMat:DimMismatch', 'Inner dimensions must agree.');
            v = kronMat.reduce(obj.terms,x,'trans',obj.trans);
        end

        function v = applyInv(obj,x)
            assert(obj.size(2) == size(x,1), 'kronMat:DimMismatch', 'Inner dimensions must agree.');
            v = kronMat.reduce(obj.terms,x,'trans',obj.trans,'inv',true);
        end

        function v = applyT(obj,x)
            assert(obj.size(2) == size(x,1), 'kronMat:DimMismatch', 'Inner dimensions must agree.');
            v = kronMat.reduce(obj.terms,x,'trans',~obj.trans);
        end

        function v = applyTinv(obj,x)
            assert(obj.size(2) == size(x,1), 'kronMat:DimMismatch', 'Inner dimensions must agree.');
            v = kronMat.reduce(obj.terms,x,'trans',~obj.trans,'inv',true);
        end

    end

    methods(Static)
        function v = reduce(terms,v,opt)
            arguments
                terms (1,:) cell
                v     {mustBeNumeric}
                opt.trans (1,1) logical = false
                opt.inv (1,1) logical = false
            end
            vSize = size(v);
            nTerms = numel(terms);
            if opt.trans
                inDims  = cellfun(@(A) size(A,1),terms);
                outDims = cellfun(@(A) size(A,2),terms);
                transFlag = 'transpose';
            else
                inDims  = cellfun(@(A) size(A,2),terms);
                outDims = cellfun(@(A) size(A,1),terms);
                transFlag = 'none';
            end
            v = reshape(v, [flip(inDims),vSize(2:end)]);
            for iTerm = nTerms:-1:1
                if opt.inv
                    v = pagemldivide(terms{iTerm},transFlag,v);
                else
                    v = pagemtimes(terms{iTerm},transFlag,v,'none');
                end
                v = permute(v,[2:nTerms,1,(nTerms+1):ndims(v)]);
            end
            v = reshape(v,[prod(outDims),vSize(2:end)]);
        end
        % TODO use tApplyOp, and optimize order
    end

end


function mustBeAllMatrices(cellArray)
for i = 1:numel(cellArray)
    assert(isnumeric(cellArray{i}) && ismatrix(cellArray{i}), ...
        'kronMat:InvalidInput','Each element in "terms" must be a numeric matrix.');
end
end