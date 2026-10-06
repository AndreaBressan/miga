classdef bSpline < aBasis
    % bSpline represents the set of univariate B-spline basis
    % corresponding to a given degree and knot vector.
    % The knot vector is not necessarily open and the basis can be
    % periodic.

    % 'degree' is the degree of the basis
    %
    % 'knots' contains the knot vector

    properties(SetAccess=protected,GetAccess=public)
        degree    (1,1) double {isinteger,mustBeNonnegative}
    end
    properties(Dependent=true)
        knots     (:,1) double
    end
    properties(Access=protected)
        compKnots (:,1) double
        exclude   (2,1) double {isinteger}
    end

    methods
        function obj = bSpline(degree,knots,periodic,domain)
            %BSPLINE Construct a B-spline basis for piecewise polynomial functions.
            %
            %   OBJ = bSpline(DEGREE,KNOTS,PERIODIC,DOMAIN)
            %   constructs a B-spline basis of polynomial degree DEGREE associated
            %   with the knot sequence KNOTS.
            %
            %   The entries of KNOTS must be in nondecreasing order. The first and
            %   last knots define the domain in the non-periodic case, while the
            %   domain can be specified independently for periodic splines.
            %
            %   At an interior knot repeated R times, consecutive polynomial pieces
            %   have DEGREE-R continuous derivatives. Knot repetitions are therefore
            %   used to control the continuity of the spline basis. The endpoints are
            %   treated as junctions with the zero function in the non-periodic case.
            %
            %   For periodic splines, the first and last knots are identified. Their
            %   repetitions are counted together, with one occurrence removed from the
            %   total count.
            %
            %   The constructor builds an augmented knot vector in OBJ.compKnots,
            %   allowing the standard de Boor algorithm to be applied on all knot
            %   spans.
            %
            %   Input arguments:
            %       DEGREE   - Polynomial degree of the B-spline basis. Must be a
            %                  nonnegative integer. Default: 0.
            %
            %       KNOTS    - Column vector containing the knots in nondecreasing
            %                  order. A knot can be repeated at most DEGREE+1 times.
            %                  In the non-periodic case, KNOTS(1) and KNOTS(end)
            %                  define the domain. Default: [-inf;inf].
            %
            %       PERIODIC - Logical flag specifying whether the B-spline basis is
            %                  periodic. Default: false.
            %
            %       DOMAIN   - Two-element vector specifying the computational domain.
            %                  For non-periodic splines, or when DOMAIN contains NaN,
            %                  the domain is set to [KNOTS(1),KNOTS(end)]. For
            %                  periodic splines, DOMAIN may be specified independently
            %                  of the first and last knots. Default: [NaN,NaN].
            %
            %   Output arguments:
            %       OBJ      - B-spline basis object containing the polynomial degree,
            %                  domain, periodicity information, basis dimension, and
            %                  the augmented knot vector required for evaluation.
            %
            %   Notes:
            %       In the non-periodic case, the knot vector is completed at both
            %       endpoints to obtain the standard open knot vector. The number of
            %       added endpoint knots is stored in OBJ.exclude.
            %
            %       In the periodic case, the knot vector is extended outside the
            %       computational domain so that the basis can be evaluated over the
            %       whole periodic domain. The resulting augmented knot vector is
            %       stored in OBJ.compKnots.
            %
            %       A periodic basis with at most one knot is reduced to a constant
            %       basis.

            arguments
                degree   (1,1) double  {isinteger} =0
                knots    (:,1) double  {}=[-inf,inf]
                periodic (1,1) logical {}=false
                domain   (1,2) double  {}=[nan,nan]
            end
            if any(isnan(domain))||~periodic; domain = [knots(1),knots(end)]; end
            simplifyToConst=periodic && numel(knots)<=1;
            obj.domain=domain;
            obj.periodic = periodic;
            obj.degree = degree*(~simplifyToConst);
            % Construct compknots, i.e. an augmented knot vector that
            % that allows to use the standard deBoor algorithm on all knot
            % spans.

            r = sum(knots==knots',2);
            assert(all(r<=obj.degree+1),"Illegal knot vector since a knot is repeated more than degree+1 times");
            r=r([1,end]);

            if ~obj.periodic
                % complete to an open knot vector and keep track of
                % added functions
                obj.exclude=degree+1-r;
                obj.compKnots = [knots(1)*ones(obj.exclude(1),1); knots; knots(end)*ones(obj.exclude(2),1)];
                obj.dimension=length(knots)-obj.degree-1;
            else
                % construct a knot vector that allows for computation on
                % the whole domain
                if knots(1)<obj.domain(1) || knots(end)>obj.domain(2)
                    error("Knots given outside of periodic domain");
                end
                reps=sum(knots==obj.domain,1);
                mrep=max(reps);
                knots=[obj.domain(1)*ones(mrep-reps(1),1);knots(knots<obj.domain(2))];
                RBefore=mod(obj.degree+1-mrep,numel(knots));
                RAfter=mod(obj.degree+1,numel(knots));
                NBefore=floor((obj.degree+1-mrep)/numel(knots))+(RBefore~=0);
                NAfter=floor((obj.degree+1)/numel(knots))+(RAfter~=0);
                st=(RBefore~=0)*numel(knots)+1-RBefore;
                en=numel(knots)*(NBefore+NAfter+1-(RAfter~=0))+RAfter;
                tk= [ diff(obj.domain).*(-NBefore:-1)+knots,knots,obj.domain(2)+(0:NAfter)*diff(obj.domain)+(knots-obj.domain(1))];
                obj.compKnots=tk(st:en);
                obj.dimension=max(1,numel(knots));
                obj.exclude=[0,0];
            end
        end

        function knots = get.knots(obj)
            %KNOTS Return the knot sequence of the B-spline basis.
            %
            %   KNOTS = obj.knots
            %   returns the knot sequence associated with the B-spline basis OBJ.
            %
            %   For a non-periodic B-spline basis, the knots added to the endpoints
            %   when constructing the augmented knot vector are removed, and only the
            %   original knot sequence is returned.
            %
            %   For a periodic B-spline basis, the augmented knot vector is restricted
            %   to the computational domain. The right endpoint is excluded so that
            %   the periodic endpoint is not represented twice.
            %
            %   Output arguments:
            %       KNOTS - Column vector containing the knot sequence associated with
            %               OBJ.


            if obj.periodic
                knots=obj.compKnots(obj.compKnots>=obj.domain(1) & obj.compKnots<obj.domain(2));
            else
                knots = obj.compKnots(1+obj.exclude(1):end-obj.exclude(2));
            end
        end

        function [sups]  = supports(obj)
            %SUPPORTS Return the supports of the B-spline basis functions.
            %
            %   SUPS = obj.supports()
            %   returns the support interval of each B-spline basis function.
            %
            %   The support of each basis function is determined by the corresponding
            %   interval in the augmented knot vector OBJ.compKnots. For a B-spline of
            %   degree P, the support of a basis function spans P+1 consecutive knot
            %   intervals.
            %
            %   For periodic B-splines, the support may extend outside the selected
            %   computational domain. In this case, the right endpoint of the support
            %   is truncated at the right endpoint of OBJ.domain.
            %
            %   Output arguments:
            %       SUPS - Array of size [1,2,OBJ.dimension]. For each basis function I,
            %              SUPS(1,1,I) and SUPS(1,2,I) contain respectively the left
            %              and right endpoints of its support.


            sups = nan(1,2,obj.dimension);
            range = (1:obj.dimension)+obj.exclude(1);
            sups(1,1,1:obj.dimension) = obj.compKnots(range);
            sups(1,2,1:obj.dimension) = min(obj.compKnots(range+obj.degree+1),obj.domain(2));
        end

        function g = greville(obj)
            %GREVILLE Return the Greville abscissae of the B-spline basis functions.
            %
            %   G = obj.greville()
            %   returns the Greville abscissae associated with the B-spline basis OBJ.
            %   For a B-spline basis of degree P, the Greville abscissa of each basis
            %   function is the average of the P knots defining its corresponding
            %   Greville point.
            %
            %   For non-periodic B-splines, the Greville abscissae are computed from
            %   the original knot sequence.
            %
            %   For periodic B-splines, the Greville abscissae are computed from the
            %   augmented knot vector. Points falling below the left endpoint of the
            %   computational domain are shifted by one period so that they are
            %   represented within the selected domain.
            %
            %   Output arguments:
            %       G - Vector containing the Greville abscissa of each B-spline basis
            %           function. Its length is equal to obj.dimension.
            %
            %   Notes:
            %       For periodic splines, the implementation assumes that the local
            %       supports do not wrap around the periodic domain more than once.


            stencil=repmat(1/obj.degree,1,obj.degree);
            if ~obj.periodic
                g=conv(obj.knots(1+(1:obj.dimension+obj.degree-1)),stencil,'valid');
            else
                g=conv(obj.compKnots(1+(1:obj.dimension+obj.degree-1)),stencil,'valid');
                inv=g<obj.domain(1);
                g(inv)=g(inv)+diff(obj.domain);
                % we are assuming that there are local functions so that supports do not wrap around the domain multiple tiems
            end
        end

        function [refBasis, refMatrix] = refine(obj,newKnots,newDeg)
            %REFINE Refine a B-spline basis by knot insertion and/or degree elevation.
            %
            %   REFBASIS = obj.refine(NEWKNOTS)
            %   constructs a refined B-spline basis by adding the knots specified in
            %   NEWKNOTS while keeping the polynomial degree unchanged. This
            %   corresponds to h-refinement.
            %
            %   REFBASIS = obj.refine([],NEWDEG)
            %   constructs a refined B-spline basis with polynomial degree NEWDEG and
            %   unchanged knot sequence. This corresponds to p-refinement.
            %
            %   REFBASIS = obj.refine(NEWKNOTS,NEWDEG)
            %   constructs the coarsest spline space containing the space spanned by
            %   OBJ and the additional knots in NEWKNOTS, with polynomial degree
            %   NEWDEG.
            %
            %   [REFBASIS,REFMATRIX] = obj.refine(...)
            %   also returns the sparse change-of-basis matrix REFMATRIX such that the
            %   coefficients of a function represented in OBJ can be transformed into
            %   coefficients with respect to REFBASIS.
            %
            %   If NEWKNOTS contains the element midpoints and NEWDEG is increased
            %   accordingly, the operation corresponds to dyadic k-refinement.
            %
            %   Input arguments:
            %       NEWKNOTS - Column vector containing the knots to be added to the
            %                  basis. For non-periodic splines, all knots must lie
            %                  within the domain. For periodic splines, knots outside
            %                  the domain are mapped periodically into the domain.
            %                  Default: [].
            %
            %       NEWDEG   - Polynomial degree of the refined basis. NEWDEG must be
            %                  an integer greater than or equal to obj.degree.
            %                  Default: obj.degree.
            %
            %   Output arguments:
            %       REFBASIS  - Refined B-spline basis of degree NEWDEG containing the
            %                   original spline space and the additional knot
            %                   refinement.
            %
            %       REFMATRIX - Sparse change-of-basis matrix mapping coefficients
            %                   from the original basis OBJ to the refined basis
            %                   REFBASIS.
            %
            %   Notes:
            %       If both NEWKNOTS is empty and NEWDEG is equal to obj.degree, the
            %       original basis is unchanged and REFMATRIX is the identity matrix.
            %
            %       For periodic B-splines, knot insertion is performed consistently
            %       across the periodic boundary so that the resulting basis remains
            %       periodic.
            %
            %       When only one output argument is requested, the change-of-basis
            %       matrix is not computed.


            arguments
                obj    bSpline
                newKnots  (:,1) double = []
                newDeg    double {mustBeNonnegative,mustBeInteger}= obj.degree
            end

            outNK=newKnots<obj.domain(1)|newKnots>obj.domain(2);
            if obj.periodic
                newKnots(outNK)= mod(newKnots(outNK)-obj.domain(1),diff(obj.domain))+obj.domain(1);
            elseif any(outNK)
                error('There are new knots are outside the domain');
            end
            degElev=newDeg-obj.degree;
            knsIns = numel(newKnots);
            if degElev<0
                error('The new degree must be greater than or equal to the spline degree');
            end
            if nargout==0;return;end
            % Computation of the new basis
            if degElev
                degKnots=repelem(unique(obj.knots),degElev);
            else
                degKnots=[];
            end
            refKnots=sort([degKnots(:);newKnots;obj.knots]);
            refBasis=bSpline(newDeg,refKnots,obj.periodic,obj.domain);
            if nargout==1;return;end

            if ~degElev && ~knsIns ; refMatrix=speye(obj.dimension); return; end
            % to compute the change of basis matrix we reduce all cases
            % to that of an open knot vector and keep track of the parts
            % to keep and discard

            tmpK=obj.compKnots;
            curMul=[find(tmpK>tmpK(1),1)-1,numel(tmpK)-find(tmpK<tmpK(end),1,'last')];
            addMul=obj.degree+1-curMul;
            openK=[repmat(tmpK(1),addMul(1),1);...
                tmpK;repmat(tmpK(end),addMul(2),1)];
            openKD=openK;
            openKR=openK;
            ignore=addMul+obj.exclude(:)';
            if degElev>0
                coefs = eye(numel(openK)-obj.degree-1);
                [degMatrix,openKD] = bSplineDegElev(obj.degree,coefs,openK,newDeg-obj.degree);
                openKD=openKD';
                refMatrix=sparse(degMatrix');
                openKR=openKD;
            else
                refMatrix=speye(numel(openK)-newDeg-1);
            end
            if knsIns>0
                if obj.periodic
                    % make sure we act correctly on basis wrapping around
                    % the periodic boundary
                    newKnots=[newKnots; newKnots+diff(obj.domain); newKnots-diff(obj.domain)];
                    newKnots=newKnots(newKnots>openKD(1)& newKnots<openKD(end));
                    newKnots=sort(newKnots);
                end
                openKR=sort([openKD;newKnots]);
                lastK=numel(openKR)-newDeg-1;
                mu=findSpans(openKD,openKR(1:lastK))-newDeg;
                c=(0:newDeg)'+mu;
                lks=openKD((1:newDeg*2)'+mu);
                v=oslo(lks,openKR((0:newDeg)'+(1:lastK)));
                r=repmat(1:lastK,newDeg+1,1);
                insMatrix = sparse(r,c,v,lastK,lastK-numel(newKnots));
                refMatrix = insMatrix*refMatrix;
            end
            refMatrix = refMatrix(ignore(1)+1:end-ignore(2),ignore(1)+1:end-ignore(2));
            if obj.periodic
                refMatrix=refMatrix(:,1:obj.dimension);
                s1=find(openKR>obj.domain,1)-newDeg-2-ignore(1);
                refMatrix(s1+refBasis.dimension +(-s1+1:0),:)= refMatrix(s1+refBasis.dimension+(-s1+1:0),:)+refMatrix(1:s1,:);
                refMatrix=refMatrix(s1+(1:refBasis.dimension),:);
            end
        end
    end


    methods (Access=protected)
        function pts=plotPoints(obj,pltdom,nPts)
            b=unique(obj.knots);
            b=b(b>=pltdom(1) & b<=pltdom(2));
            nb=numel(b);
            r=numel(b)/nPts;
            if r>=1
                pts={b(round(linspace(1,nb,nPts)))};
            else
                pts={sort([linspace(pltdom(1),pltdom(2),nPts-nb),b(:)']);};
            end
        end
        function pt=plotDecorations(obj,pltdom)
            b=unique(obj.knots);
            b=b(b>=pltdom(1) & b<=pltdom(2));
            nb=numel(b);
            pt=plot(b,zeros(size(1,nb)),'or', 'MarkerFaceColor', 'b','HandleVisibility','off');
            obj.plotDecorations@aBasis(pltdom);
        end
    end
    methods
        function [v,a] = evaluateI(obj,pts,ders)
            % returns values and actives at each point
            % one point = one column
            if iscell(pts);pts=pts{1}; end
            ord=obj.degree+1;
            nAct=min(ord,obj.dimension);
            nPts=npts(pts);
            nDer=numel(ders);
            M=numel(obj.compKnots)-obj.degree-1;
            mu=max(min(findSpans(obj.compKnots,pts(:)'),M)-obj.degree,1);
            aF=(0:obj.degree)'-obj.exclude(1)+mu;
            aF=reshape(aF,[ord,nPts])';
            if obj.periodic
                aF=1+mod(aF-1,obj.dimension);
                p=repmat((1:nPts)',nAct,1);
            else
                z=~(aF<1 | aF> obj.dimension);
                aF=aF(z);
                p=repmat((1:nPts)',ord,1);
                p=p(z);
            end
            a=gSparsity(sparse(p(:),aF(:),true(numel(p),1),nPts,obj.dimension));
            pos=a.nzPos(p(:),aF(:));

            v = cell(1,nDer);
            if nDer>0 && obj.dimension>0
                ids=(1:2*obj.degree)';
                lks=obj.compKnots(ids+mu);
                for k=numel(ders):-1:1
                    tmp=deBoor(ders(k),reshape(pts,1,[]),lks)';
                    if ~obj.periodic
                        tmp=tmp(z);
                    else
                        tmp(:,1:(obj.degree+1-obj.dimension)) = tmp(:,1:(obj.degree+1-obj.dimension))+tmp(:,obj.dimension+1:end);
                        tmp=tmp(:,1:nAct);
                    end
                    tmp(pos)=tmp;
                    v{k}=reshape(tmp,1,1,[]);
                end
            else
                v(:)={zeros(1,0,0)};
            end
        end
    end
end
