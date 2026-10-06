classdef tBasis < aBasis
    properties(GetAccess=public,SetAccess=protected)
        components      (:,1) aBasis
        vectorDimension (1,:) double
    end

    methods
        function obj = tBasis(component)
            %TBASIS Construct a tensor-product basis from univariate basis components.
            %
            %   OBJ = tBasis(COMPONENTS)
            %   constructs a tensor-product basis whose components are the basis
            %   objects specified by COMPONENTS(1,1),...,COMPONENTS(N,1).
            %
            %   The resulting basis represents the tensor product of the component
            %   spaces. Its dimension is the product of the dimensions of the
            %   component bases.
            %
            %   Input arguments:
            %       COMPONENTS - Vector of basis objects defining the components
            %             of the tensor-product basis. Each component must be
            %             univariate and must not itself be a tBasis. At most one
            %             component can be non-scalar-valued.
            %
            %   Output arguments:
            %       OBJ - Tensor-product basis containing the specified components.
            %             The main properties of OBJ are constructed from those of its
            %             components:
            %
            %             domain          - Concatenation of the component domains.
            %             vectorDimension - Dimensions of the component bases.
            %             dimension       - Product of the component dimensions.
            %             periodic        - Periodicity of each component.
            %             rangeSize       - Concatenation of the component range sizes.
            %
            %   Notes:
            %       All components must be univariate, i.e., their domain must have
            %       one row.
            %
            %       Nested tensor-product bases are not supported: none of the
            %       components can be an instance of tBasis.
            %
            %       At most one component may be non-scalar-valued. This restriction
            %       prevents multiple components from having a non-scalar range.


            arguments (Repeating)
                component (:,1) aBasis {}
            end
            obj.components=cat(1,component{:});
            obj.domain=cat(1,obj.components.domain);
            obj.vectorDimension=[obj.components.dimension];
            obj.dimension=prod(obj.vectorDimension);
            obj.periodic=cat(1,obj.components.periodic);
            obj.rangeSize=cat(2,obj.components.rangeSize);
            assert(sum(arrayfun(@(c)prod(c.rangeSize),obj.components)>1)<=1,...
                "At most one component of a tBasis can be non-scalar-valued");
            assert(all(arrayfun(@(c)~isa(c,'tBasis'),obj.components)),'Nested tensor bases is not supported.');
            % TODO remove the following limitation so that one could do FEM
            % in space IGA in time or similar hybrids
            assert(all(arrayfun(@(c)size(c.domain,1),obj.components)==1),'All components must be univariate.');
        end

        function [w,a]= evaluateCombinationI(obj, coefs, pts, ders)
            if issparse(coefs) || ~iscell(pts)
                [w,a]=evaluateCombinationI@aBasis(obj,coefs,pts,ders);
                return;
            end

            isGrid=iscell(pts);
            dim=size(obj.domain,1);
            if ~isGrid
                nPts=npts(pts);
                pts=mat2cell(pts,ones(size(obj.domain,1),1));
            else
                uNpts=cellfun(@npts,pts);
                nPts=prod(uNpts);
            end

            nDer=size(ders(:,:),2);
            nFun=size(coefs,2);
            nRan=prod(obj.rangeSize);
            uNfun=obj.vectorDimension;
            if obj.dimension==0; w=nan(nRan,nDer,nPts,nFun); return; end

            [nT,rMap]=tNormalizeTable(ders);
            cG=tReductionNetwork(nT);

            [uv,ua]=arrayfun(@(c)obj.components(c).evaluateI(pts{c},rMap(c).operators),...
                1:dim,'UniformOutput',false);
            mat=@(c,o)full( ua{c}.toSparse(uv{c}{o}))';
            for c=dim:-1:1
                nOp=numel(rMap(c).operators);
                mats=arrayfun(mat,repmat(c,1,nOp),1:nOp,'UniformOutput',false);
                uOps{c}=cellfun(@(m) @(t) tApplyOp(t,m,c+1), mats,'UniformOutput',false );
            end
            arg=reshape(coefs,[1,uNfun,nFun,nRan]);
            a=fSparsity([nPts,nFun]);
            F=cG.FconstructO2M(uOps);
            w=cellfun(@(m)reshape(m,[nRan,1,nPts*nFun]),F(arg),'UniformOutput',false);
        end

        function [b]  = supports (obj)
            %SUPPORTS Return the supports of the tensor-product basis functions.
            %
            %   B = obj.supports()
            %   returns the support of each basis function of the tensor-product basis
            %   OBJ.
            %
            %   The support of a tensor-product basis function is obtained by combining
            %   the supports of its component basis functions. In particular, for a
            %   tensor-product basis with D components, the support is represented by
            %   D intervals, one for each component direction.
            %
            %   Output arguments:
            %       B - Array of size [D,2,OBJ.vectorDimension], where D is the
            %           dimension of the domain. For each tensor-product basis
            %           function, B(:,1,:) and B(:,2,:) contain respectively the lower
            %           and upper bounds of its support in each component direction.
            %           The last dimensions are flattened so that each combination of
            %           component basis functions corresponds to one tensor-product
            %           basis function.
            %
            %   Notes:
            %       The supports of the tensor-product basis functions are obtained
            %       directly from the supports of the component bases. The ordering of
            %       the tensor-product basis functions follows the ordering induced by
            %       the component dimensions.


            dim=size(obj.domain,1);
            b=zeros([dim,2,obj.vectorDimension]);
            ind=cat(2,{0,[1,2]},repmat({':'},1,dim));
            sz=[1,2,ones(1,dim)];
            for c=1:numel(obj.components)
                ind{1}=c;
                sz(c+2)=obj.vectorDimension(c);
                sc=obj.components(c).supports;
                b(ind{:})=b(ind{:})+reshape(sc,sz);
                sz(c+2)=1;
            end
            b=reshape(b,dim,2,[]);
        end

        function S = sparsity(obj,other,sdom)
            %SPARSITY Construct the sparsity pattern between two basis objects.
            %
            %   S = obj.sparsity()
            %   returns the sparsity pattern describing the overlap between the basis
            %   functions of OBJ.
            %
            %   S = obj.sparsity(OTHER)
            %   returns the sparsity pattern between the basis functions of OBJ and
            %   OTHER.
            %
            %   S = obj.sparsity(OTHER,SDOM)
            %   restricts the overlap computation to the subdomain specified by SDOM.
            %   SDOM contains one interval for each spatial dimension.
            %
            %   When OBJ and OTHER are both tensor-product bases, the sparsity pattern
            %   is constructed from the sparsity patterns of the corresponding
            %   component bases. In particular, the sparsity of each component is
            %   computed independently on the corresponding component of SDOM, and
            %   the resulting component sparsity patterns are combined into a
            %   tensor-product sparsity pattern.
            %
            %   Input arguments:
            %       OTHER - Basis object whose basis functions are compared with those
            %               of OBJ. If omitted, OTHER is set equal to OBJ. When OTHER
            %               is a tBasis, its components are matched with the
            %               corresponding components of OBJ.
            %
            %       SDOM  - Subdomain on which the overlap between basis functions is
            %               considered. SDOM must have size [D,2], where D is the
            %               dimension of the domain. Each row specifies the lower and
            %               upper bounds of the subdomain in the corresponding
            %               direction. Default: [-inf,inf] in every direction.
            %
            %   Output arguments:
            %       S - aSparsity object describing the nonzero interactions between
            %           the basis functions of OBJ and OTHER. If both bases are
            %           tensor-product bases, S is a tSparsity object obtained from
            %           the component sparsity patterns.
            %
            %   Notes:
            %       If OTHER is not a tBasis, the sparsity computation is delegated to
            %       the implementation in aBasis.
            %
            %       For tensor-product bases, the component bases are paired according
            %       to their ordering, and the corresponding components of SDOM are
            %       used for the individual sparsity computations.


            arguments
                obj
                other (1,1) aBasis=obj
                sdom  (:,2) double=repmat([-inf,inf],size(obj.domain,1),1)
            end
            if isa(other,'tBasis')
                comps=arrayfun(@(i)obj.components(i).sparsity(other.components(i),sdom(i,:)),1:size(obj.domain,1));
                S=tSparsity(comps);
            else
                S=obj.sparsity@aBasis(other,sdom);
            end
        end
    end
    methods(Access=protected)
        function pts=plotPoints(obj,pltdom,nPts)
            cpts=@(c)obj.components(c).plotPoints(pltdom(c,:),nPts(c,:));
            pts=arrayfun(cpts,1:numel(obj.components));
        end
    end
    methods
        function [v,a]=evaluateI(obj,pts,ders)
            dim=size(obj.domain,1);
            ran_c=arrayfun(@(b)prod(b.rangeSize),obj.components);
            ran=prod(ran_c);

            % Compute points
            isGrid=iscell(pts);
            if isGrid
                nPts_c=cellfun(@npts,pts);
                nPts=prod(nPts_c);
            else
                nPts=npts(pts);
                nPts_c=nPts*ones(1,dim);
                pts=mat2cell(pts,ones(size(obj.domain,1),1));
            end

            % Decide what to compute
            [nT,rMap]=tNormalizeTable(ders);
            [v_c,a_c]=arrayfun(@(c)obj.components(c).evaluateI(pts{c},rMap(c).operators),1:dim,'UniformOutput',false);
            cG=tReductionNetwork(nT);

            if isGrid
                % TODO change the requirement of tSparsity.toSparse so that
                % arbitrary component sparsity are accepted and the
                % transformation to col major is not needed

                % convert univariate results to col major
                tmp_conv1=arrayfun(@(c)a_c{c}.toSparse((1:ran_c(c))'+ran_c(c)*(0:a_c{c}.nnz-1)),1:dim,"UniformOutput",false );
                ord_c=cellfun(@nonzeros,tmp_conv1,"UniformOutput",false);
                A_c=cellfun(@(mat)gSparsity(mat~=0),tmp_conv1,"UniformOutput",false);
                V_c=arrayfun(@(c) cellfun(@(data) data(:,1,ord_c{c}),v_c{c},"UniformOutput",false),1:dim,"UniformOutput",false);

                % Use the converted sparsities
                a=tSparsity([A_c{:}]);
                mat=@(c,o) reshape(V_c{c}{o},[ran_c(c),ones(1,c-1),A_c{c}.nnz,ones(1,dim-c)]);
                for c=dim:-1:1
                    nOp=numel(rMap(c).operators);
                    mats=arrayfun(mat,repmat(c,1,nOp),1:nOp,'UniformOutput',false);
                    uOps{c}=cellfun(@(m) @(t) t.*m, mats,'UniformOutput',false );
                end
                F=cG.FconstructO2M(uOps);
                v=cellfun(@(m)reshape(m,[ran,1,a.nnz]),F(1),'UniformOutput',false);
            else
                % A similar strategy cannot work here because the non-zeros
                % are a different number for each point and we do not have
                % a Kronecker structure of the evaluation matrices as for grid

                % In this case the sparsity is the concatenation of pointwise
                % outer product and flattening.
                % To exploit the tensor product structure the algorithm
                % classifies points based on number of active

                nnz_pt_c=cellfun(@(a)a.nnzRow,a_c,'UniformOutput',false);
                nnz_pt_pc=cat(2,nnz_pt_c{:});
                [nnz_tp,~,nnz_tp_ind]=unique(nnz_pt_pc,'rows','first');
                dim_mul=cumprod([1,obj.vectorDimension]);
                v=arrayfun(@(o)zeros(ran,1,0),1:size(ders,2),"UniformOutput",false);
                act=zeros(0,1);
                pts=act;
                for tp=1:size(nnz_tp,1)
                    cur_pts= find(nnz_tp_ind==tp);
                    cur_nPts= numel(cur_pts);
                    cur_nAct_c=nnz_tp(tp,:);
                    cur_nAct=prod(cur_nAct_c);
                    cur_act=1;
                    for c=dim:-1:1
                        sp_mat=a_c{c}.toSparse(1:a_c{c}.nnz);
                        [~,cur_act_c,cur_ord]=find(sp_mat(cur_pts,:)); % this is column major
                        mat=@(dat) reshape(dat(:,1,cur_ord),[ran_c(c),cur_nPts,ones(1,c-1),cur_nAct_c(c)]);
                        mats=cellfun(mat,v_c{c},'UniformOutput',false);
                        uOps{c}=cellfun(@(m) @(t) t.*m, mats,'UniformOutput',false );
                        tmp_act=reshape(cur_act_c,[cur_nPts, ones(1,c-1), cur_nAct_c(c)]);
                        cur_act=cur_act+(tmp_act-1)*dim_mul(c);
                    end
                    %cur_act=permute(cur_act,[2:dim+1,1]);
                    F=cG.FconstructO2M(uOps);
                    cur_v=cellfun(@(m)reshape(m,[ran,1,cur_nAct*numel(cur_pts)]),F(1),'UniformOutput',false);
                    v=cellfun(@(a,b)cat(3,a,b),v,cur_v',"UniformOutput",false);
                    pts=cat(1,pts,repmat(cur_pts(:),cur_nAct,1));
                    act=cat(1,act(:),cur_act(:));
                end
                a=gSparsity(sparse(pts,act,1:numel(pts),nPts,obj.dimension));
            end

        end
    end
end
