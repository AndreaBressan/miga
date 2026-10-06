function geo_b = boundary(index, geo) 
%BOUNDARY Restrict a geometry to one of its boundary faces.
%
%   GEO_B = boundary(INDEX,GEO)
%   constructs the geometry induced on the boundary face specified by
%   INDEX of the geometry GEO.
%
%   The boundary index follows the convention
%
%       2*d-1  -> lower boundary in direction d
%       2*d    -> upper boundary in direction d
%
%   The returned geometry has one less parameter dimension than GEO, while
%   retaining the same range dimension.
%
%   Input arguments:
%       INDEX - Positive integer identifying the boundary face. If DIR is
%               the associated parameter direction, then
%
%                   DIR = ceil(INDEX/2).
%
%               Odd indices select the lower boundary and even indices
%               select the upper boundary.
%
%       GEO   - Geometry to be restricted to the selected boundary.
%               GEO must provide the fields domain, domainDim, rangeDim,
%               evaluate, and jacobian.
%
%   Output arguments:
%       GEO_B - Geometry restricted to the selected boundary. Its
%               properties are:
%
%               domainDim - GEO.domainDim-1.
%               rangeDim  - GEO.rangeDim.
%               domain    - GEO.domain with the coordinate corresponding
%                           to the selected boundary removed.
%               evaluate  - Function handle for evaluating the restricted
%                           geometry.
%               jacobian  - Function handle for evaluating the restricted
%                           Jacobian.


dir=ceil(index/2);
ate= 2-mod(index,2);
sdir=setdiff(1:geo.domainDim,dir);

geo_b.domainDim=geo.domainDim-1;
geo_b.rangeDim=geo.rangeDim;

geo_b.domain=geo.domain(sdir,:);
geo_b.evaluate=@(varargin)evaluateGeo(geo,dir,sdir,ate,varargin);%geo.evaluate(insert(geo.domain(dir,ate),dir,x));
geo_b.jacobian=@(x)restrictJ(geo, dir,sdir, ate, x);

end

function varargout=evaluateGeo(geo,dir,sdir,ate,args)
varargout = cell(max(numel(args)-1,1),1);
if numel(args)==1
    args{2}='value';
end
for i=2:numel(args)
    switch class(args{i})
        case 'char'
            switch args{i}
                case 'value'
                    varargout{i-1}=geo.evaluate(insert(geo.domain(dir,ate),dir,args{1}));
                case 'jacobian'
                    varargout{i-1}=restrictJ(geo, dir,sdir, ate, args{1});
            end
    end
end
end



function y=insert(val,i,x)
if iscell(x)
    y=cat(1,x{1:i-1},{val},x{i:end,:});
else
    sx=size(x); 
    sy=sx+[1,zeros(numel(sx)-1)];
    p=prod(sx(2:end));
    y=reshape(cat(1,x(1:i-1,:),repmat(val,1,p),x(i:end,:)),sy);
end
end

function J=restrictJ(geo,dir,sdir,ate,x)
T=geo.jacobian(insert(geo.domain(dir,ate),dir,x));
sT=size(T);
J=reshape(T(:,sdir,:), [sT(1), sT(2)-1, sT(3:end)]);
end
