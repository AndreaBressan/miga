classdef exprBasis < aBasis
    % exprBasis represents a function or a basis that is described by
    % a string where x1,...,xdim are the coordinates dimensions

    properties(GetAccess=public,SetAccess=protected)
        expression       string {}
        symbolicFun      sym    {}
        matlabFun        struct {}
        variablesSym     sym    {}
        variablesNames   string {}
        parametersNames  string {}
        parametersValue  double {}
    end

    methods
        function obj = exprBasis(domain,expressions,rangeSize,periodic, options)
            arguments
                domain (:,2) double {} = [-inf,inf]
                expressions string = "1"
                rangeSize (:,1) double  = ones(0,1)
                periodic  (:,1) logical = []
                options.parameters (1,2) cell ={string.empty,double.empty}
            end
            obj.rangeSize=rangeSize;
            obj.domain=domain;
            dim=size(obj.domain,1);
            if ~isempty(periodic); obj.periodic=periodic; else, obj.periodic=false(dim,1);end
            obj.dimension=numel(expressions)/prod(rangeSize);
            obj.expression=reshape(expressions,[],1);
            obj.symbolicFun=str2sym(obj.expression);
            obj.variablesSym=sym('x',[dim,1]);
            obj.variablesNames=string(obj.variablesSym);
            obj.parametersNames=string(setdiff(symvar(obj.symbolicFun),obj.variablesSym));
            obj.parametersValue=nan(size(obj.parametersNames));
            obj=obj.setParameters(options.parameters{1},options.parameters{2});
        end

        function [b] = supports(obj)
            b = repmat(obj.domain,1,1,obj.dimension);
        end

        function obj=setParameters(obj,params, values)
            arguments
                obj    (1,1) exprBasis
                params (:,1) string = string.empty
                values (:,1) double = double.empty
            end
            for i=1:numel(params)
                pos=find(strcmp(obj.parametersNames,params(i)));
                if ~isempty(pos)
                    obj.parametersValue(pos)=values(i);
                else
                    warning("Parameter '%s' is not in the expression", params(i));
                end
            end
            if any(isnan(obj.parametersValue))
                obj.matlabFun=struct("expression",{} ,"function", {});
            else
                expr=obj.symbolicFun;
                for i=1:numel(obj.parametersNames)
                    expr=simplify(subs(expr,obj.parametersNames(i),obj.parametersValue(i)));
                end
                cfh=functionFromExpression(expr,obj.variablesNames);
                rangeDim=prod(obj.rangeSize);
                functionDim=obj.dimension;
                fh=@(varargin)exprBasis.reshapeResult(cfh,rangeDim,functionDim,varargin);
                obj.matlabFun=struct("expression", {expr}, "function",{fh});
            end
        end
    end

    methods
        function [v,a] = evaluateI(obj,pts,ders)
            dim=size(obj.domain,1);
            nPts=npts(pts);
            nw=size(ders(:,:),2);
            if iscell(pts)
                ptsl=cellfun(@numel,pts);
                GF=@(i)reshape(...
                    repmat( pts{i}(:)',prod(ptsl(1:i-1)),1,prod(ptsl(i+1:dim)))...
                    , [1,1,nPts]);
                xs=arrayfun(GF,1:dim,"UniformOutput",false);
            else
                xs=num2cell(pts(:,:),2);
            end
            v=cell(1,nw);
            for w=1:size(ders(:,:),2)
                curPD=ders(:,w);
                curFH=obj.handle(curPD);
                cv=curFH(xs{:});
                v{w}=reshape(cv,[prod(obj.rangeSize),1,nPts*obj.dimension]);
            end
            a=fSparsity([nPts,obj.dimension]);
        end

        function fh = handle(obj,pd)
            dim=size(obj.domain,1);
            fsz=size(obj.matlabFun,1:dim);
            pdc=num2cell(pd+1);
            if all(fsz(:)>=pd+1) && ~isempty(obj.matlabFun(pdc{:}).expression)
                fh=obj.matlabFun(pdc{:}).function;
            else
                fh=obj.computeMatlabFunction(pd);
            end
        end


        function fh = computeMatlabFunction(obj,pd)
            dim=size(obj.domain,1);
            % resize cache array and fill required item
            fsz=size(obj.matlabFun,1:max(dim,2));
            if isscalar(pd); pdf=[pd,0]; else, pdf=pd'; end
            if any(fsz<=pd) % regrow cell if necessary
                fsz=max(fsz,pdf+1);
                mpos=num2cell(fsz);
                obj.matlabFun(mpos{:})=struct("expression",{string.empty},"function",{function_handle.empty});
            end
            mpos=num2cell(pd+1);
            pos=sub2ind(fsz,mpos{:});

            derVar=num2cell(obj.variablesSym(repelem(1:dim,pd)));
            expr=diff(obj.matlabFun(1).expression,derVar{:});
            cfh=functionFromExpression(expr,obj.variablesNames);
            rangeDim=prod(obj.rangeSize);
            functionDim=obj.dimension;
            fh=@(varargin)exprBasis.reshapeResult(cfh,rangeDim,functionDim,varargin);
            obj.matlabFun(pos).expression=expr;
            obj.matlabFun(pos).function=fh;
        end

    end
    methods(Static)
        function out=reshapeResult(cfh,rangeDim,functionDim,args)
            comps=cellfun(@(f) f(args{:}),reshape(cfh,rangeDim,1,1,functionDim),'UniformOutput',false);
            out=cell2mat(comps);
        end
    end
end
