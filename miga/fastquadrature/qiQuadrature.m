classdef qiQuadrature < aIntegration

    properties
        map
        QI
        reference
    end

    methods
        function obj=qiQuadrature(QI,refQ,domain)
            arguments
                QI         (:,:) struct
                refQ       (1,1) aIntegration
                domain     (:,2) double = zeros(0,2)
            end
            if isempty(domain)
                obj.domain=QI.basis.domain;
            else
                obj.domain=domain;
            end
            obj.nodes=QI.nodes;
            obj.reference=refQ;
            obj.QI=QI;
            obj.map=obj.QI.basis.evaluate(refQ.nodes)*QI.mat;
            
            % compute weights (value,domain,rangeSize,dimension)
            refBasisInt=obj.reference.integrate(constantFunction(1,obj.domain,[],1),obj.QI.basis);
            obj.weights=refBasisInt'*QI.mat;
        end

        function [iOp]=getTestOperator  (obj,test,request)
            nReq=size(request,2);
            iOp=cell(nReq,1);
            for r=1:nReq
                refOp=obj.reference.integrate(constantFunction(1,obj.domain,[],1),...
                    test,cat(3,request(:,r),zeros(size(request,1),1)),obj.QI.basis);
                iOp{r}=(refOp*obj.QI.mat)';
            end
        end

        function [iOp]=getMatrixOperator(obj,test,trial,request,sparsity)
            nReq=size(request,2);
            iOp=cell(nReq,1);
            for r=1:nReq
                refOp=obj.reference.getMatrixOperator(test,trial,request(:,r,:),sparsity);
                iOp{r}=obj.map'*refOp{1};
            end
        end


        function [iOp]=getApplyOperator(obj,test,trial,request)
            %TODO make something more sensible for QI, this is basically a weighted quadrature
            nReq=size(request,2);
            iOp=cell(nReq,1);
            for r=1:nReq
                ref=obj.reference.integrate(constantFunction(1,obj.domain,[],1), test,request(:,r,:),trial);
                val=trial.evaluate(obj.nodes, parDers(request(:,r,2)));                
                iOp{r}=lsqminnorm(val',ref');
            end
        end

    end
end