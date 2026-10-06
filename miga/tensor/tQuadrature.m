classdef tQuadrature < aIntegration    
    properties
        components
    end
    
    methods
        function obj = tQuadrature(components)
            %TQUADRATURE Construct a tensor-product quadrature rule.
            %
            %   OBJ = tQuadrature(COMPONENTS)
            %   constructs a tensor-product quadrature rule from the quadrature rules
            %   specified in COMPONENTS.
            %
            %   Each component defines the quadrature nodes and weights in one
            %   coordinate direction. The resulting quadrature combines the component
            %   rules into a multidimensional tensor-product quadrature.
            %
            %   Input arguments:
            %       COMPONENTS - Column vector of aIntegration objects defining the
            %                    component quadrature rules. Each component corresponds
            %                    to one direction of the tensor-product quadrature.
            %                    Default: aIntegration.empty.
            %
            %   Output arguments:
            %       OBJ - Tensor-product quadrature object. Its main properties are
            %             constructed from the component quadratures:
            %
            %             domain     - Concatenation of the component domains.
            %             nodes      - Tensor-product quadrature nodes.
            %             weights    - Tensor-product quadrature weights.
            %             components - Component quadrature objects.
            %
            %   Notes:
            %       The nodes and weights of the component quadratures are converted
            %       to cell arrays when necessary and combined to form the
            %       multidimensional tensor-product representation.


            arguments
                components (:,1) aIntegration = aIntegration.empty
            end
            obj.domain=cat(1,components.domain);
            tmp1=arrayfun( @(c)ensureCell(c.nodes), components,'UniformOutput',false);
            obj.nodes=[tmp1{:}];
            tmp2=arrayfun( @(c)ensureCell(c.weights), components,'UniformOutput',false);
            obj.weights=[tmp2{:}];
            obj.components=components;
        end
        
        function [iOp] = getIntegrateOperator(obj)
            dirOp=arrayfun(@(c)c.getIntegrateOperator(),obj.components,'UniformOutput',false);
            dirNN=cellfun(@npts,{obj.components.nodes});
            function out=tOp(in)
                sz=size(in,1:max(3,ndims(in)));
                in=reshape(in,[sz(1),dirNN,sz(3:end)]);
                out=in;
                for dir=1:numel(obj.components)
                    out=tApplyOp(out, dirOp{dir},dir+1);
                end
                out=reshape(out,sz([1,3:numel(sz)]));
            end
            iOp=@tOp;
        end

        function [iOp] = getTestOperator(obj,test,request)
            if isa(test,'tBasis')
                [req,map]=tNormalizeTable(request);
                optimized=tReductionNetwork(req);
                uOps=arrayfun(@(q,uB,um)q.getTestOperator(uB,um.operators),...
                        obj.components,test.components,map,'UniformOutput',false);
                dim=numel(obj.components);
                nTerm=size(request,2);
                startingDim=cellfun(@npts,{obj.components.nodes});
                for c=1:dim
                    uOps{c}=cellfun(@(m)@(t)tApplyOp(t,m,c+1),uOps{c},'UniformOutput',false );
                end
                integrator=optimized.FconstructM2O(uOps);
                iOp=@(m)reshape(integrator(prepare(m,startingDim,nTerm)),test.dimension,[]);
            else
                iOp=obj.getTestOperator@aIntegration(test,request);
            end
        end

        function [iOp]=getMatrixOperator(obj,test,trial,request,sparsity)
            if isa(test,'tBasis') && isa(trial,'tBasis') && isa(sparsity, "tSparsity")
                [req,map]=tNormalizeTable(request);
                optimized=tReductionNetwork(req);
                uOps=arrayfun(@(q,ute,utr,um,sp)q.getMatrixOperator(ute,utr,um.operators,sp),...
                        obj.components,test.components,trial.components,map,sparsity.uSparsity,'UniformOutput',false);
                dim=numel(obj.components);
                nTerm=size(request,2);
                startingDim=cellfun(@npts,{obj.components.nodes});
                for c=1:dim
                    uOps{c}=cellfun(@(m) @(t)tApplyOp(t,m,c+1), uOps{c},'UniformOutput',false );
                end
                integrator=optimized.FconstructM2O(uOps);
                trf=prod(test.rangeSize)*prod(trial.rangeSize);
                iOp=@(m)reshape(integrator(prepare(m,startingDim,nTerm)),size(m,1)*trf,sparsity.nnz,[]);
            else
                iOp=obj.getMatrixOperator@aIntegration(test,trial,request,sparsity);
            end
        end

        function [iOp]=getApplyOperator(obj,test,trial,request)
            if isa(test,'tBasis') && isa(trial,'tBasis')
                [req,map]=tNormalizeTable(request);
                optimized=tReductionNetwork(req);
                uOps=arrayfun(@(q,ute,utr,um,sp)q.getApplyOperator(ute,utr,um.operators),...
                        obj.components,test.components,trial.components,map,'UniformOutput',false);
                dim=numel(obj.components);
                nTerm=size(request,2);
                startingDim=cellfun(@npts,{obj.components.nodes});
                for c=1:dim
                    uOps{c}=cellfun(@(m) @(t)tApplyOp(t,m,c+1), uOps{c},'UniformOutput',false );
                end
                integrator=optimized.FconstructM2O(uOps);
                trf=prod(test.rangeSize)*prod(trial.rangeSize);
                iOp=@(m) reshape(integrator(prepare(m,startingDim,nTerm)),size(m,1)*trf,test.dimension,[]);
            else
                iOp=obj.getApplyOperator@aIntegration(test,trial,request);
            end
        end

    end
end

function cell=prepare(m,startingDim,nTerm)
m2ca=cat(2,num2cell([size(m,1),startingDim]),{ones(nTerm,1)});
cell=mat2cell(reshape(m,[size(m,1),startingDim,size(m,3)]),m2ca{:});
end


function out=ensureCell(arg)
% TODO maybe this is not necessary
% ensureCell, if the input arg is a cellarray leave it unmodified otherwise
% wraps it in a cell array
if iscell(arg)
    out=arg;
else
    out={arg};
end
end
