function s=surfLines(varargin)
% surfsmooth wraps standard surf function and by default uses 
% interpolation shader, EdgeColor=None, and adds high resolution grid
% lines to the surfaces as separate objects.
narginchk(1,inf)
linePos=find(strcmp(varargin,'GridLineNumber'),1);
edgePos=find(strcmp(varargin,'EdgeColor'),1);
facePos=find(strcmp(varargin,'FaceColor'),1);
parentPos=find(strcmp(varargin,'Parent'),1);
displayNamePos=find(strcmp(varargin,'DisplayName'),1);

if isempty(facePos);   varargin=cat(2,varargin,{'FaceColor','interp'});end
if isempty(edgePos)
    hasEdge=false;
    varargin=cat(2,varargin,{'EdgeColor','none'});
elseif strcmp(varargin{edgePos+1},'none')
    hasEdge=false;
else
    hasEdge=true;
end

if isempty(linePos)
    % default to 0 lines if edges are drown
    lineNum=(~hasEdge)*[11,11]; 
else
    lineNum=varargin{linePos+1};
    varargin([linePos,linePos+1])=[];
end
if isscalar(lineNum); lineNum =lineNum*[1,1]; end

if ~isempty(parentPos)
    hg = hggroup(varargin{parentPos+1});
    varargin{parentPos+1}=hg;
else
    hg = hggroup(gca);
    varargin=cat(2,varargin,{'Parent',hg});
end
hg.Annotation.LegendInformation.IconDisplayStyle='on';
hg.HandleVisibility='on';
if isempty(displayNamePos)
    hg.DisplayName='my data';
else
    hg.DisplayName=varargin{displayNamePos+1};
end

s=surf(varargin{:});
X=s.XData;
Y=s.YData;
Z=s.ZData;
sz=size(X);
xpos=round(linspace(1,sz(2), lineNum(1)));
ypos=round(linspace(1,sz(1), lineNum(2)));
holdState=ishold;
hold on
for i =xpos; plot3(X(:,i),Y(:,i),Z(:,i),'-k','Parent',hg); end
for i =ypos; plot3(X(i,:),Y(i,:),Z(i,:),'-k','Parent',hg); end
if~holdState; hold off; end
end