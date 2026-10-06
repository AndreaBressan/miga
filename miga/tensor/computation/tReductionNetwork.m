classdef tReductionNetwork
    % tReductionNetwork represents a set of tensor contractions
    % that either start from a tensor and output multiple tensors,
    % such as the evaluation of multiple partial derivatives
    % of a tensor basis, or given in input multiple tensors the
    % output is
    properties(GetAccess=public,SetAccess=protected)
        compTable
    end

    properties(Access=protected)
        % resorted table
        sCompTab
        sTabOrd
        sTabOrdInv
        % algorithm internals
        nTmp (1,1) double
        nNod (1,1) double
        % node description
        nodeComp
        nodeOp
        nodeTerm
        % operation edges
        eExt1
        eExt2
    end

    methods
        function obj = tReductionNetwork(compTable)
            %TREDUCTIONNETWORK Construct an instance of this class
            %   Detailed explanation goes here
            arguments
                compTable (:,:) double {}=[]
            end
            obj.compTable=compTable;
            dim=size(compTable,1);
            nTerm=size(compTable,2);

            [obj.sCompTab,obj.sTabOrdInv]=sortrows([zeros(nTerm,1),compTable',(1:nTerm)']);
            delim=(obj.sCompTab(2:end,:)-obj.sCompTab(1:end-1,:))~=0;
            for c=2:dim+2; delim(:,c)=delim(:,c)|delim(:,c-1); end
            nTmp=sum(delim,1)+1;
            obj.nNod=sum(nTmp);
            nIdsD=cumsum([ones(1,dim+2);delim],1)+cumsum([0,nTmp(1:end-1)]);
            [~,nodeTerm]=unique(nIdsD);
            obj.nodeTerm=mod(nodeTerm-1,nTerm)+1;

            opCel=arrayfun(@(c)obj.sCompTab([true;delim(:,c)],c),1:dim+2,'UniformOutput',false);
            obj.nodeOp=cat(1,opCel{:});
            obj.nodeComp=repelem(0:dim+1,nTmp)';
            obj.sTabOrd(obj.sTabOrdInv)=1:nTerm;
            obj.eExt1=(2:obj.nNod)';
            eExCe=arrayfun(@(c)nIdsD([true;delim(:,c)],c-1),2:dim+2,'UniformOutput',false);
            obj.eExt2=cat(1,eExCe{:});
        end

        function G=graphO2M(obj)
            G=obj.graph('O2M');
        end

        function G=graphM2O(obj)
            G=obj.graph('M2O');
        end

        function G=graph(obj,type)
            dim=size(obj.compTable,1);
            function fstr=fStrNode(c)
                if c==0
                    fstr=string(['N_{%d} ( *',repmat(', *',1,dim-c-1),')']);
                elseif c<=dim
                    fstr=string(['N_{%d} (',repmat('% 2d',c>0,1),repmat(',% 2d',1,c-1),repmat(', *',1,dim-c),')']);
                else
                    fstr=string(['N_{%d} (% 2d',repmat(',% 2d',1,dim-1),') [%d]']);
                end
            end
            nodeLabel=arrayfun(@(c,id)string(sprintf(fStrNode(c),[id,obj.sCompTab(obj.nodeTerm(id),2:c+1)])),obj.nodeComp,(1:obj.nNod)');
            edgeLabel=arrayfun(@(id,c,op)string(sprintf('E_{%d} (dir %d, op %d)',id,c,op)),(1:obj.nNod-1)', obj.nodeOp(2:end),obj.nodeComp(2:end));
            nodeTable=table(nodeLabel, obj.nodeComp, obj.nodeOp, 'VariableNames',{'Name' 'Component' 'Operator'});
            switch type
                case 'O2M'
                    edgeTable=table([obj.eExt2,obj.eExt1],edgeLabel, obj.nodeOp(2:end), obj.nodeComp(2:end), 'VariableNames',{'EndNodes'  'EdgeLabel' 'Operator' 'Component'});
                case 'M2O'
                    edgeTable=table([obj.eExt1,obj.eExt2],edgeLabel, obj.nodeOp(2:end), obj.nodeComp(2:end), 'VariableNames',{'EndNodes'  'EdgeLabel' 'Operator' 'Component'});
                otherwise, error('Unknown type of reduction');
            end
            G=digraph(edgeTable,nodeTable);
        end

        function [FO2M]=FconstructO2M(obj,ops)
            c  = obj.nodeComp;
            op = obj.nodeOp;
            dim=size(obj.compTable,1);
            CO2M=cell(obj.nNod,1);
            for j=obj.nNod:-1:1
                o_nodes=obj.eExt1(obj.eExt2==j);
                if c(j)>dim % out
                    CO2M{j}=@(in){in};
                elseif c(j)>=1
                    CO2M{j}=@(in) executeAll(ops{c(j)}{op(j)}(in),CO2M(o_nodes));
                else
                    CO2M{j}=@(in) executeAll(in,CO2M(o_nodes));
                end
            end
            % reorder output
            tabOrd=obj.sTabOrd;
            FO2M=@(in) subsref( CO2M{1}(in),struct('type',{'()'},'subs',{{tabOrd}}));
        end

        function [FM2O]=FconstructM2O(obj,ops)
            c  = obj.nodeComp;
            op = obj.nodeOp;
            dim=size(obj.compTable,1);
            CM2O=cell(obj.nNod,1);
            for j=obj.nNod:-1:1
                i_nodes=obj.eExt1(obj.eExt2==j);
                if c(j)>dim % in
                    assert(numel(i_nodes)==0);
                    CM2O{j}=@(in)in{op(j)};
                elseif c(j)>0
                    sumF=sumAll(CM2O(i_nodes));
                    CM2O{j}=@(in) ops{c(j)}{op(j)}(sumF(in));
                else
                    CM2O{j}=sumAll(CM2O(i_nodes));
                end
            end
            FM2O=CM2O{1};
        end

    end
end

function [out]=executeAll(arg,funs)
tmp=cellfun(@(f)f(arg),funs(:),'UniformOutput',false);
out=cat(1,tmp{:});
end

function [out]=sumAll(funs)
nel=numel(funs);
tmp=funs;
while nel>1
    for j=1:floor(nel/2)
        tmp(j)={@(in)tmp{2*j-1}(in)+tmp{2*j}(in) };
    end
    if mod(nel,2)==1
        j=j+1;
        tmp(j)=tmp(nel);
    end
    nel=j;
end
out=tmp{1};
end