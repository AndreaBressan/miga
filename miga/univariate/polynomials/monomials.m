classdef monomials <aBasis
    % monomials represents the polynomial basis given by monomials

    % 'degree' contains the degree of the basis

    properties
        degree
    end

    methods
        function obj = monomials(degree,domain,periodic)
            arguments
                degree   (1,1) double =-1;
                domain   (1,2) double =[-inf,inf]
                periodic (1,1) logical=false
            end
            obj.degree=degree;
            obj.dimension=obj.degree+1;
            obj.domain=domain;
            obj.periodic=periodic;
            obj.rangeSize=[];
        end
        function [b]  = supports     (obj)
            b=repmat(obj.domain,[1,1,obj.degree+1]);
        end
    end
    methods
        function [v,a]=evaluateI(obj,pts,ders)
            if iscell(pts); pts=pts{1}; end
            pts = pts(:);
            nPts = numel(pts);
            nAct = obj.degree + 1;
            nDer = numel(ders);
            
            a = fSparsity([nPts, nAct]);
            v = cell(1, nDer);
            
            maxDer=max(ders);
            curDer=pts.^(0:obj.degree);
            curPos=1;
            curDeg=1:obj.degree;
            for w = 0:maxDer
                if ders(curPos)==w
                    v{curPos}=reshape(curDer,1,1,nAct*nPts);
                    curPos=curPos+1;
                end
                curDer(:,2:end)=curDer(:,1:end-1).*curDeg;
                curDer(:,1)=0;
                curDeg=max(curDeg-1,0);
            end
        end
    end
end