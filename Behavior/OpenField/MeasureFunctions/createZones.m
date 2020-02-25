function [zones] = createZones(origin, diameter,varargin)
%This function creates pie shpaed quadrants for beahvioral testing in circular apparatus
%by LB & MG March 2019

%Input:
%       origin: 1x2 matrix of xy coords for origin, ex: [0,0] for
%       coordinates centered at 0,0.
%       diameter: scalar quantity that reflects the diameter of the behavioral apparatus.
%Output:
%       zones: a nx2 matrix of xy coordinates for quadrant verticies.
%       Column 1 = x, column 2 = y. These can be used as an input to
%       computeDwell for calculating dwell time.
%Vargains:
%       type: general shape of quadrant either pie or annulus, pie is
%       default.
%       numQuad: number of symmetrical quadrants, default is 4.
%       annulusSize: Size of annulus based on the percentage of the
%       apparatus circumference.
%       fig: plots a figure of

p = inputParser;

% %ADD ERROR MESSAGE AND FILTER
%%errorMsg = 'Value must be either type or pie as a string argument';
% % validationFcn_char = @(x) assert(~ischar(x),errorMsg);

% addParameter(p,'type','pie',validationFcn_char);

addParameter(p,'type','pie');
addParameter(p,'numQuad',4);
addParameter(p,'annulusSize',.8);
addParameter(p,'fig',0);
addParameter(p,'title','Quadrants');

% type=parse(p,Parameters,'');
p.parse(varargin{:});

type=p.Results.type;
numQuad=p.Results.numQuad;
annulusSize=p.Results.annulusSize;
fig=p.Results.fig;
Figtitle=p.Results.title;

%Params for quadrants
x0=origin(:,1);%%origin
y0=origin(:,2);%%origin
r=diameter/2;%%radius
n=numQuad; %number of quadrants
tempX=(cos(linspace(-pi,pi,1000))*r)+x0;
tempY=(sin(linspace(-pi,pi,1000))*r)+y0;

if strcmp(type,'pie')
    
    zones(:,1)=(r+r*.25)*cos(linspace(-pi,pi,n+1))+x0; %sets zones %25 outside of enviornmental boundary to account for error
    zones(:,2)=(r+r*.25)*sin(linspace(-pi,pi,n+1))+y0;
    
else
    
    zones(:,1)=tempX*annulusSize;
    zones(:,2)=tempY*annulusSize;
    
end

if fig==1 && strcmp(type,'pie')
    
    f=figure;
    f.Color=[1 1 1];
    tet=linspace(-pi,pi,n+1); %%allows for even space between each line??
    xi=r*cos(tet)+x0;
    yi=r*sin(tet)+y0;
    title(Figtitle)
    for k=1:numel(xi)
        plot([x0 xi(k)],[y0 yi(k)],'k','LineWidth',2); hold on
    end
    
    hold on;
    
    plot(tempX,tempY,'k','LineWidth',2);
    
    axis image
    axis off
    
elseif fig==1 && strcmp(type,'annulus')
    f=figure;
    f.Color=[1 1 1];
    title(Figtitle)
    
    plot(tempX,tempY,'-k'); hold on; plot(zones(:,1),zones(:,2));
    
    axis image
    axis off
    
end

end
