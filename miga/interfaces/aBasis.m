classdef aBasis < matlab.mixin.Heterogeneous
    % aBasis is the interface that allows to evaluate a set of functions:
    % a single function, a basis or generating set.

    properties(GetAccess=public,SetAccess=protected)
        % 'domain' and 'periodic' describe the domain of the function[s].
        %
        % 'domain' is a n x 2 matrix describing a n-dimensional
        %          axis aligned cuboid in R^n. The i-th row represents the
        %          i-th axis, while the columns are respectively the vertex
        %          having minimal and maximal coordinate.
        %
        % 'periodic' is a n x 1 logical vector describing if the domain is
        %          periodi in each dimension. This affects how the value are
        %          computed at points outside the domain.
        %          0: non-periodic direction
        %          1: periodic direction
        %
        % 'dimension' is the number of functions modelled by the basis
        %
        % 'rangeSize' describes the size of the range.
        %             Scalars are represented by []
        %             [m] are m-dimensional vectors
        %             [a1,...,am] are order m tensors with sizes a1, ..., am

        domain    (:,2) double
        periodic  (:,1) logical
        dimension (1,1) double
        rangeSize (1,:) double
    end

    properties
        % outOfDomainError sets the behavior of evaluation at points outside of
        % the domain:
        % - if outOfDomainError is true a call to evaluate with a point
        %   outside of the domain will raise an exception
        % - if outOfDomainError is false[the default] a 0 value will be
        %   returned for all basis functions and all points outside the domain.
        outOfDomainError (1,1) logical = false
    end

    properties(GetAccess=protected,SetAccess=private)
        derivativeShift (:,1) double =0
    end

    methods(Abstract)
        % EVALUATEI implements the evaluation in subclasses.
        % The interface is similar to that of evaluate, but
        % - ders must be specified and it must be a list of partial derivatives
        % - v is a cell with one partial derivative per element shaped as
        %   [prod(rangeSize), numDerivatives, num nonzeros]
        [v,a] = evaluateI(obj,pts,ders)
    end
    methods(Abstract)
        % SUPPORTS returns a obj.domainDim x 2 x obj.dimension tensor with the bounding
        % boxes of the supports of the functions. If the domain is periodic then
        % the boxes can extend outside of fundamental domain and the user should
        % take care of periodicity and the identification of points.
        [b]  = supports     (obj)
    end

    methods
        function [pts,vps] = validatePoints(obj,pts)
            %VALIDATEPOINTS Validate and preprocess evaluation points.
            %
            %   [PTS,VPS] = obj.validatePoints(PTS)
            %   checks the evaluation points with respect to the domain of OBJ,
            %   removes points outside non-periodic directions, and maps points
            %   outside periodic directions back into the corresponding domain.
            %
            %   The function supports two representations of the evaluation points:
            %
            %       1. Cell representation:
            %          PTS is a cell array containing one vector for each coordinate
            %          direction. In this case, the points form a Cartesian grid.
            %
            %       2. Explicit-coordinate representation:
            %          PTS is a double array with one row for each coordinate
            %          direction. Each column represents one evaluation point.
            %
            %   For non-periodic directions, points outside the domain are discarded.
            %   If the property OBJ.outOfDomainError is true, an error is instead
            %   raised when any point lies outside the domain.
            %
            %   For periodic directions, points outside the domain are mapped back
            %   into the domain using periodicity. In particular, a coordinate is
            %   mapped to the corresponding point in the interval
            %
            %       [obj.domain(d,1), obj.domain(d,2)].
            %
            %   Input arguments:
            %       OBJ  - aBasis object defining the computational domain and its
            %              periodic directions.
            %
            %       PTS  - evaluation points, either as:
            %              * a cell array containing one vector per coordinate
            %                direction (Cartesian grid); or
            %              * a double matrix with one row per coordinate direction
            %                and one column per point.
            %
            %   Output arguments:
            %       PTS  - validated and possibly modified evaluation points.
            %              Points outside non-periodic directions are removed, while
            %              points outside periodic directions are mapped back into
            %              the domain.
            %
            %       VPS  - logical vector indicating which original evaluation points
            %              belong to the domain in all non-periodic directions.
            %              For a Cartesian grid, VPS refers to the points obtained
            %              from the Cartesian product of the coordinate vectors.
            per=obj.periodic;
            chk=~per;
            dirs=1:size(obj.domain,1);
            function nv=setInDomain(v,d)
                nv=obj.domain(d,1)+mod(v-obj.domain(d,1),obj.domain(d,2)-obj.domain(d,1));
            end
            if iscell(pts) % cartesian grid
                assert(numel(pts)==size(obj.domain,1),'Invalid points, maybe wrong dimension?')
                vps=true;
                for d=dirs(chk)
                    ok=pts{d}>=obj.domain(d,1) & pts{d}<= obj.domain(d,2);
                    assert(~obj.outOfDomainError || all(ok),'Evaluation outside of the domain requested');
                    vps=vps & reshape(ok,[ones(1,d-1) numel(ok) 1]);
                    pts{d}=pts{d}(ok);
                end
                nRep=ones(1,numel(dirs));
                nRep(obj.periodic)=cellfun(@npts,pts(obj.periodic));
                vps=repmat(vps,nRep);
                vps=vps(:);
                for d=dirs(per)
                    ok=pts{d}>=obj.domain(d,1) & pts{d}<= obj.domain(d,2);
                    pts{d}(~ok)=setInDomain(pts{d}(~ok),d);
                end
            else % explicit coordinates
                assert(size(pts,1)==size(obj.domain,1),'MIGA:POINTS_WROND_DIMENSION','Invalid points, maybe wrong dimension?')
                vps= pts(:,:)>=obj.domain(:,1) & pts(:,:)<=obj.domain(:,2);
                assert(~obj.outOfDomainError || all(vps,'all'),'MIGA:POINTS_OUTSIDE_DOMAIN','Evaluation outside of the domain requested');
                for d=dirs(per)
                    ok=vps(d,:);
                    pts(d,~ok)=setInDomain(pts(d,~ok),d);
                end
                vps=all(vps(chk,:),1);
                pts=pts(:,vps);
            end
        end

        function [pds,func] = parseEvaluationRequest(obj,request)
            %PARSEEVALUATIONREQUEST Parse and organize an evaluation request.
            %
            %   [PDS,FUNC] = obj.parseEvaluationRequest(REQUEST)
            %   parses an evaluation request and converts it into a set of partial
            %   derivative requests PDS together with a function FUNC that reconstructs
            %   the requested output quantities from the evaluated derivatives.
            %
            %   REQUEST can be specified in one of the following forms:
            %
            %       1. An evaluationOp object:
            %          The derivative information required by the operator is
            %          obtained by calling its WHAT method. FUNC is then the
            %          evaluation method associated with the operator.
            %
            %       2. A cell array containing evaluationOp objects and/or explicit
            %          derivative requests. Each entry is parsed recursively.
            %          Duplicate derivative requests are identified and evaluated
            %          only once. FUNC then combines the computed quantities to
            %          produce the outputs corresponding to the original requests.
            %
            %   Input arguments:
            %       OBJ     - aBasis object. Its domain dimension is used to interpret
            %                 the derivative requests.
            %
            %       REQUEST - evaluation request. It can be either an evaluationOp
            %                 object or a cell array containing valid evaluation
            %                 requests.
            %
            %   Output arguments:
            %       PDS     - matrix containing the partial derivative multi-indices
            %                 required to evaluate REQUEST. Each column represents
            %                 one partial derivative request.
            %
            %       FUNC    - function handle used to construct the requested output
            %                 quantities from the evaluated derivatives. The function
            %                 has the form
            %
            %                     FUNC(VC,DOMAINDIM,RANGESIZE)
            %
            %                 where VC contains the evaluated derivatives,
            %                 DOMAINDIM is the domain dimension, and RANGESIZE is the
            %                 range size of the basis.
            domainDim=size(obj.domain,1);

            if isa(request,'evaluationOp')
                pds=request.what(domainDim);
                func=@request.evaluate;
            elseif isa(request,'cell')
                [s_pds,s_func]=cellfun(@(s_req)obj.parseEvaluationRequest(s_req),request,'UniformOutput',false);
                blkSize       =cellfun(@(ds)size(ds,2),s_pds);
                [pdst,~,dI]=unique(cat(2,s_pds{:})','rows');
                pds_ind=mat2cell(dI,blkSize);
                pds    =pdst';
                func   = @(vc,domainDim,rangeSize) cellfun(@(f,w)f(vc(w),domainDim,rangeSize),s_func(:),pds_ind(:),'UniformOutput',false);
            else
                error('Invalid evaluation request: they can be explicit partial derivatives or named operators from evaluationOps folder or cells containing any of the two');

            end
        end
    end

    methods

        function [varargout] = evaluate(obj,pts,ders)
            %EVALUATE Evaluate basis functions and their derivatives at given points.
            %
            %   [V] = obj.evaluate(PTS,DERS)
            %   evaluates the basis functions associated with OBJ at the points PTS
            %   and returns the quantities specified by the derivative requests
            %   DERS.
            %
            %   [V,A] = obj.evaluate(PTS,DERS)
            %   additionally returns A, which contains information on the active
            %   basis functions at the evaluation points. This output is available
            %   when one more output argument than the number of derivative requests
            %   is requested.
            %
            %   Input arguments:
            %       OBJ   - aBasis object defining the basis to be evaluated.
            %
            %       PTS   - evaluation points. The points must satisfy the
            %               requirements imposed by the basis domain and by the
            %               mustBeValidPts validation function.
            %
            %       DERS  - derivative/evaluation requests (vector). Each request specifies
            %               the derivative of the basis functions to be evaluated.
            %               If omitted, the basis functions themselves are evaluated.
            %
            %   Output arguments:
            %        V     - values of the requested basis functions or derivatives
            %               evaluated at PTS. The output format depends on the
            %               requested derivatives and on the basis dimension/range.
            %
            %       A     - active-basis information associated with the evaluation
            %               points. This output is returned only when explicitly
            %               requested as an additional output argument.
            arguments
                obj aBasis
                pts {mustBeValidPts(obj, pts)}
            end
            arguments (Repeating)
                ders
            end

            domainDim=size(obj.domain,1);
            nRequest=numel(ders);
            nOut=max(1,nRequest);
            nPts=npts(pts);
            returnActive= (nargout==nOut+1);

            if ~returnActive && nargout>0 && nargout~=nOut; error("Requested %d outputs from the arguments, but %d assigned", nOut, nargout); end
            if isempty(ders); ders = {parDers(zeros(size(obj.domain,1),1))}; end

            [pds,func]  =obj.parseEvaluationRequest(ders);
            [vpts,valid]=obj.validatePoints(pts);
            [v,a] =obj.evaluateI(vpts,pds+obj.derivativeShift);
            if any(~valid); a=rSparsity(a,find(valid),[nPts,obj.dimension]); end

            if ~returnActive
                for i=1:numel(v)
                    v{i}=a.toFull(v{i},3);
                end
            end

            outs=func(v,domainDim,obj.rangeSize);
            if returnActive
                varargout=cat(1,outs(:),{a});
            else
                varargout=outs(:);
            end
        end


        function [varargout] = evaluateCombination(obj,coefs,pts,ders)
            %EVALUATECOMBINATION Evaluate a linear combination of basis functions.
            %
            %   [V] = obj.evaluateCombination(COEFS,PTS,DERS)
            %   evaluates the linear combination of the basis functions defined by
            %   the coefficients COEFS at the points PTS. The quantities specified
            %   by the derivative requests DERS are returned.
            %
            %   [V,A] = obj.evaluateCombination(COEFS,PTS,DERS)
            %   additionally returns A, which contains information on the active
            %   basis functions at the evaluation points.
            %
            %   Input arguments:
            %       OBJ   - aBasis object defining the basis to be evaluated.
            %
            %       COEFS - coefficients defining the linear combination of the
            %               basis functions. The size and structure of COEFS must be
            %               compatible with the basis represented by OBJ.
            %
            %       PTS   - evaluation points. The points must satisfy the
            %               requirements imposed by the basis domain and by the
            %               mustBeValidPts validation function.
            %
            %       DERS  - derivative/evaluation requests. Each request specifies
            %               the derivative of the linear combination to be evaluated.
            %               If omitted, the linear combination itself is evaluated.
            %
            %   Output arguments:
            %       V  - values of the requested quantities evaluated at PTS.
            %               The output format depends on the requested derivatives
            %               and on the dimension and range of the basis.
            %
            %       A     - active-basis information associated with the evaluation
            %               points. This output is returned only when explicitly
            %               requested as an additional output argument.
            arguments
                obj   aBasis
                coefs double
                pts {mustBeValidPts(obj, pts)}
            end
            arguments (Repeating)
                ders (:,1) 
            end
            domainDim=size(obj.domain,1);
            numRequest=numel(ders);
            numOut=max(1,numRequest);
            returnActive= (nargout==numOut+1);

            if ~returnActive && nargout~=numOut; error("Requested %d outputs from the arguments, but %d assigned", numOut, nargout); end
            if isempty(ders); ders = {parDers(zeros(size(obj.domain,1),1))}; end

            [pds,func]  =obj.parseEvaluationRequest(ders);
            [vpts,valid]=obj.validatePoints(pts);
            [v,a]       =obj.evaluateCombinationI(coefs,vpts,pds);
            if any(~valid); a=rSparsity(a,find(valid),[nPts,obj.dimension]); end
 
            outs=func(v,domainDim,obj.rangeSize);
            if returnActive
                varargout=cat(1,outs(:),{a});
            else
                varargout=outs(:);
            end
        end

        function [w,ac] = evaluateCombinationI(obj,coefs,pts,ders)
            %EVALUATECOMBINATIONI Evaluate a linear combination of basis functions.
            %
            %   [W,AC] = obj.evaluateCombinationI(COEFS,PTS,DERS)
            %   evaluates the linear combination of the basis functions defined by
            %   OBJ, with coefficients COEFS, at the points PTS and for the
            %   derivatives specified by DERS.
            %
            %   This is the internal evaluation routine used by EVALUATECOMBINATION.
            %   It operates directly on the internal representation of the basis
            %   functions and exploits their sparsity whenever possible.
            %
            %   Input arguments:
            %       OBJ   - aBasis object defining the basis functions.
            %
            %       COEFS - matrix containing the coefficients of the linear
            %               combination. Each row corresponds to one basis function,
            %               hence
            %
            %                   size(COEFS,1) = obj.dimension.
            %
            %               The columns of COEFS identify the different linear
            %               combinations to be evaluated.
            %
            %       PTS   - evaluation points.
            %
            %       DERS  - derivative requests. Each column of DERS specifies a
            %               derivative of the basis functions to be evaluated.
            %
            %   Output arguments:
            %       W     - cell array containing the requested evaluations. Each
            %               cell corresponds to one derivative request. The entries
            %               are arranged according to the active basis functions and
            %               have size
            %
            %                   [prod(obj.rangeSize), 1, AC.nnz].
            %
            %       AC    - sparsity structure describing the active linear
            %               combinations at the evaluation points. Its shape is
            %
            %                   [nPts, size(COEFS,2)].
            %
            %               The nonzero entries identify the pairs of evaluation
            %               points and coefficient columns for which the linear
            %               combination has an active contribution.
            assert(size(coefs,1)==obj.dimension,"The coefficients must be a matrix having one row per basis function");
            nPts=npts(pts);
            nDder=size(ders(:,:),2);
            [v,a]= obj.evaluateI(pts,ders);
            basisFull=isa(a,'fSparsity');
            if issparse(coefs) && ~basisFull
                ac=gSparsity(a.toSparse(true(a.nnz,1),1)*coefs~=0);
                [rs,cs]=ac.pos;
                extract=sub2ind(ac.shape,rs,cs);
                for i=1:nDder
                    tmp=a.toSparse(v{i},3)*coefs;
                    v{i}=full(tmp(extract));
                end
            else
                coefs=full(coefs);
                ac=fSparsity([nPts,size(coefs(:,:),2)]);
                if basisFull
                    opFun=@(i)reshape(a.toFull(v{i},3),prod(obj.rangeSize)*nPts,obj.dimension)';
                else
                    opFun=@(i)a.toSparse(v{i},3)';
                end
                for i=1:nDder
                    v{i}=tApplyOp(coefs,opFun(i),1);
                end
            end
            w=cellfun(@(m)reshape(m,[prod(obj.rangeSize),1,ac.nnz]),v,"UniformOutput",false);
        end

        function [B] = derivedGenerators(obj,amount)
            %DERIVEDGENERATORS Construct a basis for derivatives of the basis functions.
            %
            %   B = obj.derivedGenerators(AMOUNT)
            %   returns a basis object B whose evaluation corresponds to derivatives
            %   of the basis functions of OBJ. This construction is useful, for
            %   example, when discretizing spaces of differential forms.
            %
            %   Input arguments:
            %       OBJ    - aBasis object defining the original basis.
            %
            %       AMOUNT - column vector specifying the derivative order in each
            %                coordinate direction. Its size must be equal to the
            %                spatial dimension of the domain:
            %
            %                   size(AMOUNT) = [dim,1],
            %
            %                where dim = size(obj.domain,1).
            %
            %                The entries of AMOUNT specify the additional derivative
            %                order applied in each coordinate direction.
            %
            %   Output arguments:
            %       B      - aBasis object representing the derived basis. The
            %                derivative shift of OBJ is increased by AMOUNT, so that
            %                subsequent evaluations of B correspond to the
            %                corresponding derivatives of the original basis
            %                functions.
            assert(iscolumn(amount) && size(amount,1)==size(obj.domain,1),"Amount should be a column of the proper size");
            assert(all(amount+obj.derivativeShift>=0),"The order of derivative should be positive")
            B=obj;
            B.derivativeShift=amount+obj.derivativeShift;
        end

        function [S] = sparsity(obj,other,sdom)
            %SPARSITY Determine the sparsity pattern induced by overlapping supports.
            %
            %   S = obj.sparsity(OTHER)
            %   returns the sparsity pattern describing which basis functions of OBJ
            %   and OTHER have overlapping supports.
            %
            %   S = obj.sparsity(OTHER,SDOM)
            %   restricts the overlap test to the subdomain SDOM.
            %
            %   The rows of S correspond to the basis functions of OBJ, while the
            %   columns correspond to the basis functions of OTHER. An entry of S is
            %   nonzero if the corresponding supports overlap in the prescribed
            %   domain.
            %
            %   Input arguments:
            %       OBJ   - basis object whose functions define the rows of S.
            %
            %       OTHER - basis object whose functions define the columns of S.
            %               By default, OTHER is set equal to OBJ.
            %
            %       SDOM  - subdomain used to restrict the overlap test. It is
            %               specified as a matrix with one row for each coordinate
            %               direction and two columns containing the lower and upper
            %               bounds:
            %
            %                   SDOM = [xmin xmax;
            %                           ymin ymax;
            %                           ...       ].
            %
            %               By default, SDOM is the whole domain, represented by
            %               [-inf,inf] in every coordinate direction.
            %
            %   Output arguments:
            %       S     - gSparsity object representing the sparsity pattern
            %
            arguments
                obj
                other  (1,1) aBasis=obj
                sdom   (:,2) double=repmat([-inf,inf],size(obj.domain,1),1)
            end
            ms=obj.supports();
            os=other.supports();
            assert(size(ms,1)==size(os,1), "Domain dimensions do not agree");
            os=reshape(os,[size(os,[1,2]),1,size(os,3)]);
            NNZ=~any(ms(:,2,:)<sdom(:,1) | ms(:,2,:)<=os(:,1,:,:) |  ms(:,1,:)>sdom(:,2) |  ms(:,1,:)>=os(:,2,:,:),1);
            NNZ=sparse(reshape(NNZ,size(NNZ,[3,4])));
            S=gSparsity(NNZ);
        end
    end

    % ================== METHODS for GRAPHICS and DEBUGGING ==============
    methods
        function ax = plot(obj,nPts,pd)
            %PLOT Plot the basis functions over the computational domain.
            %
            %   AX = obj.plot()
            %   plots all the basis functions defined by OBJ using a default number
            %   of evaluation points.
            %
            %   AX = obj.plot(NPTS)
            %   specifies the number of evaluation points used in each coordinate
            %   direction. NPTS can be either a scalar, in which case the same
            %   number of points is used in every direction, or a column vector
            %   containing one value for each coordinate direction.
            %
            %   AX = obj.plot(NPTS,PD)
            %   plots the quantities obtained by applying the evaluation operator PD
            %   to the basis functions. By default, PD corresponds to the evaluation
            %   of the basis functions themselves.
            %
            %   Input arguments:
            %       OBJ   - basis object to be plotted.
            %
            %       NPTS  - number of evaluation points in each coordinate direction.
            %               It can be either:
            %
            %                   * a scalar, used for all coordinate directions;
            %                   * a column vector with one entry per coordinate
            %                     direction.
            %
            %               The default value is 150 points in each direction.
            %
            %       PD    - evaluation operator used to evaluate the basis functions.
            %               The default is the identity evaluation operator, i.e.,
            %               the basis functions themselves are plotted.
            %
            %   Output arguments:
            %       AX    - handle to the current axes containing the plot.
            arguments
                obj
                nPts (:,1) double = 150
                pd   (:,1) evaluationOp = parDers(zeros(size(obj.domain,1),1))
            end

            if prod(obj.rangeSize)~=1; error("Cannot plot non scalar valued functions"); end
            dim=size(obj.domain,1);
            if isscalar(nPts); nPts=nPts*ones(dim,1); end
            if numel(nPts)~=dim; error("The number of points should be either a scalar or a vector of length domain dimension"); end
            extraWidth=diff(obj.domain,1,2)/6;
            plotDomain=[obj.domain(:,1)-extraWidth,obj.domain(:,2)+extraWidth];
            pts=obj.plotPoints(plotDomain,nPts);
            v = obj.evaluate(pts,pd);
            names=arrayfun(@(b)sprintf("$B_{%d}$",b),1:obj.dimension,'UniformOutput',false);
            holdState=ishold; hold on;
            switch dim
                case 1
                    p=plot(pts{1},squeeze(v));
                    set(p(:),{'DisplayName'},names(:));
                case 2
                    X=repmat(pts{1}(:),1,nPts(2));
                    Y=repmat(pts{2}(:)',nPts(1),1);
                    v=reshape(v,[nPts(:)',obj.dimension]);
                    for b=1:obj.dimension; surfLines(X,Y,v(:,:,b),'DisplayName',names{b},'FaceColor','interp'); end
                    view(3);
                otherwise
                    error("Cannot plot basis for domains of dimension >2")
            end
            obj.plotDecorations(plotDomain);
            legend('Interpreter','latex','ItemHitFcn',@cb_legend,'FontSize',20);
            if ~holdState; hold off; end
            ax=gca;
        end
    end


    methods(Access=protected)
        function pts=plotPoints(obj,pltdom,Npts)
            dim=size(obj.domain,1);
            upt=@(c)sort([linspace(pltdom(c,1),pltdom(c,2),Npts(c)-2),obj.domain(c,:)]);
            pts=arrayfun(@(c)upt(c),1:dim,"UniformOutput",false);
        end
        function pt=plotDecorations(obj,pltdomain)
            dim=size(obj.domain,1);
            switch dim
                case 1
                    pt=xline(obj.domain(obj.domain>=pltdomain(1)&obj.domain<=pltdomain(2)),'--b','HandleVisibility','off');
                case 2

                otherwise

            end
        end
    end



    methods(Static,Sealed,Access = protected)
        % necessary function to allow allocation of arrays with invalid objects
        function default=getDefaultScalarElement()
            default=aBasisDefault();
        end
    end
end


function cb_legend(~,evt)
if strcmp(evt.Peer.Visible,'on')
    evt.Peer.Visible = 'off';
else
    evt.Peer.Visible = 'on';
end
end

function mustBeValidPts(obj, pts)
if isnumeric(pts)
    % Case 1: double array with domainDim rows
    if ~isa(pts, 'double') || size(pts,1) ~= size(obj.domain,1)
        error("pts must be a double array with %d rows.", size(obj.domain,1));
    end
elseif iscell(pts)
    % Case 2: cell array of 1D double vectors (any length)
    for i = 1:numel(pts)
        pt = pts{i};
        if ~isa(pt, 'double') || ~isvector(pt)
            error("Element %d of pts must be a 1D double vector (any length).", i);
        end
    end
else
    error("pts must be a double array or a cell array of 1D double vectors.");
end
end
