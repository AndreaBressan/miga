classdef lagrange <aBasis
    % lagrange represents the basis of Lagrange polynomials
    % that are polinomials that are dual to the evaluation
    % operators at given knots

     % 'knots' contains the knot vector
     %
     % 'degree' contains the degree of the basis

    properties(GetAccess=public,SetAccess=protected)
        knots
        degree
    end

    methods
        function obj = lagrange(knots,domain,periodic)
            arguments
                knots    (:,1) double =[];
                domain   (1,2) double =[-inf,inf]
                periodic (1,1) logical=false
            end
            obj.knots=knots;
            obj.degree=numel(obj.knots)-1;
            obj.dimension=obj.degree+1;
            obj.domain=domain;
            obj.periodic=periodic;
            obj.rangeSize=[];
        end
        function [b]  = supports     (obj)
            b=repmat(obj.domain,[1,1,obj.degree+1]);
        end
    end
    methods(Access=protected)
        function pt=plotDecorations(obj,pltdomain)
            b=unique(obj.knots(obj.knots>=pltdomain(1) & obj.knots<=pltdomain(2)));
            nb=numel(b);  
            pt=plot(b,zeros(size(1,nb)),'or', 'MarkerFaceColor', 'b','HandleVisibility','off');
            obj.plotDecorations@aBasis(pltdomain);
        end
    end
    methods
        function [v,a]=evaluateI(obj,pts,ders)
            if iscell(pts);pts=pts{1}; end
            pts=pts(:)';
            ord=obj.degree+1;
            a=fSparsity([numel(pts),ord]);
            nDer=size(ders,2);
            nPts=size(pts,2);
            v=cell(1,nDer);
            
            kDiffs=obj.knots'-obj.knots;
            kDiffs(logical(eye(ord)))=1;
            normF=prod(kDiffs,2);
            pDiffs=-obj.knots+pts;

            for w=1:nDer
                der=ders(w);
                nPro=obj.degree-der;
                if nPro<0
                    v{w}=zeros(1,1,a.nnz);
                elseif nPro==0
                    % choosek(v,k) does not return an empty array when k=0
                    t= (factorial(der)*(-1)^obj.degree).*normF.\ones(ord,1,1,nPts);
                    v{w}= reshape( permute(t,[4,1,2,3]),[1,1,a.nnz]);
                else
                    c=nchoosek(1:ord-1,nPro);
                    nAdd=size(c,1);
                    c=reshape(mod((0:ord-1)'+c(:)',ord)+1,ord,nAdd,nPro);
                    tDiff=reshape(pDiffs(c,:),[ord,nAdd,nPro,nPts]);
                    t=(factorial(der)*(-1)^obj.degree).*normF.\sum(prod(tDiff,3),2);                
                    v{w}= reshape(permute(t,[4,1,2,3]),[1,1,a.nnz]) ;
                end
            end
        end
    end
end