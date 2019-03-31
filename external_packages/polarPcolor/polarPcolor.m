function [varargout] = polarPcolor(R,theta,Z,varargin)
% [h,c] = polarPcolor1(R,theta,Z,varargin) is a pseudocolor plot of matrix 
% Z for a vector radius R and a vector angle theta. 
% The elements of Z specify the color in each cell of the 
% plot. The goal is to apply pcolor function with a polar grid, which 
% provides a better visualization than a cartesian grid.
%
%% Syntax
% 
% [h,c] = polarPcolor(R,theta,Z)
% [h,c] = polarPcolor(R,theta,Z,'Ncircles',10)
% [h,c] = polarPcolor(R,theta,Z,'Nspokes',5)
% [h,c] = polarPcolor(R,theta,Z,'Nspokes',5,'colBar',0) 
% [h,c] = polarPcolor(R,theta,Z,'Nspokes',5,'labelR','r (km)')
% 
% INPUT
%	* R :
%        - type: float
%        - size: [1 x Nrr ] where Nrr = numel(R).
%        - dimension: radial distance.
%	* theta : 
%        - type: float
%        - size: [1 x Ntheta ] where Ntheta = numel(theta).
%        - dimension: azimuth or elevation angle (deg).
%        - N.B.: The zero is defined with respect to the North.
%	* Z : 
%        - type: float
%        - size: [Ntheta x Nrr]
%        - dimension: user's defined .
%	* varargin:
%        - Ncircles: number  of circles for the grid definition.
%        - Nspokes: number of spokes for the grid definition.
%        - colBar: display the colorbar or not.
%        - labelR: legend for R.
% 
% 
% OUTPUT
% h: returns a handle to a SURFACE object.
% c: returns a handle to a COLORBAR object.
%
%% Examples 
% R = linspace(3,10,100);
% theta = linspace(0,180,360);
% Z = linspace(0,10,360)'*linspace(0,10,100);
% figure
% polarPcolor(R,theta,Z,'Ncircles',3)
%
%% Author
% Etienne Cheynet, University of Stavanger, Norway. 28/05/2016
% see also pcolor
% 

%%  InputParseer
p = inputParser();
p.CaseSensitive = false;
p.addOptional('Ncircles',5);
p.addOptional('Nspokes',8);
p.addOptional('labelR','');
p.addOptional('colBar',1);
p.parse(varargin{:});

Ncircles = p.Results.Ncircles ;
Nspokes = p.Results.Nspokes ;
labelR = p.Results.labelR ;
colBar = p.Results.colBar ;
%% Preliminary checks
% case where dimension is reversed
Nrr = numel(R);
Noo = numel(theta);
if isequal(size(Z),[Noo,Nrr]),
    Z=Z';
end

% case where dimension of Z is not compatible with theta and R
if ~isequal(size(Z),[Nrr,Noo])
    fprintf('\n')
    fprintf([ 'Size of Z is : [',num2str(size(Z)),'] \n']);
    fprintf([ 'Size of R is : [',num2str(size(R)),'] \n']);
    fprintf([ 'Size of theta is : [',num2str(size(theta)),'] \n\n']);
    error(' dimension of Z does not agree with dimension of R and Theta')
end
%% data plot
rMin = min(R);
rMax = max(R);
thetaMin=min(theta);
thetaMax =max(theta);
% Definition of the mesh
Rrange = rMax - rMin; % get the range for the radius
rNorm = R/Rrange; %normalized radius [0,1]
% get hold state
cax = newplot;
% transform data in polar coordinates to Cartesian coordinates.
YY = (rNorm)'*cosd(theta);
XX = (rNorm)'*sind(theta);
% plot data on top of grid
h = pcolor(XX,YY,Z,'parent',cax);
shading flat
set(cax,'dataaspectratio',[1 1 1]);axis off;
if ~ishold(cax);
    % make a radial grid
    hold(cax,'on')
    % Draw circles and spokes
    createSpokes(thetaMin,thetaMax,Ncircles,Nspokes);
    createCircles(rMin,rMax,thetaMin,thetaMax,Ncircles,Nspokes)
end

%% PLot colorbar if specified
if colBar==1,
    c =colorbar('location','WestOutside');
    caxis([quantile(Z(:),0.01),quantile(Z(:),0.99)])
else
    c = [];
end

%% Outputs
nargoutchk(0,2)
if nargout==1,
    varargout{1}=h;
elseif nargout==2,
    varargout{1}=h;
    varargout{2}=c;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nested functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function createSpokes(thetaMin,thetaMax,Ncircles,Nspokes)
        
        circleMesh = linspace(rMin,rMax,Ncircles);
        spokeMesh = linspace(thetaMin,thetaMax,Nspokes);
        contour = abs((circleMesh - circleMesh(1))/Rrange+R(1)/Rrange);
        cost = cosd(90-spokeMesh); % the zero angle is aligned with North
        sint = sind(90-spokeMesh); % the zero angle is aligned with North
        for kk = 1:Nspokes
            plot(cost(kk)*contour,sint(kk)*contour,'k:',...
                'handlevisibility','off');
            % plot graduations of angles
            % avoid superimposition of 0 and 360
            if and(thetaMin==0,thetaMax == 360),
                if spokeMesh(kk)<360,
                    
                    text(1.05.*contour(end).*cost(kk),...
                        1.05.*contour(end).*sint(kk),...
                        [num2str(spokeMesh(kk),3),char(176)],...
                        'horiz', 'center', 'vert', 'middle');
                end
            else
                text(1.05.*contour(end).*cost(kk),...
                    1.05.*contour(end).*sint(kk),...
                    [num2str(spokeMesh(kk),3),char(176)],...
                    'horiz', 'center', 'vert', 'middle');
            end
            
        end
    end
    function createCircles(rMin,rMax,thetaMin,thetaMax,Ncircles,Nspokes)
        
        % define the grid in polar coordinates
        angleGrid = linspace(90-thetaMin,90-thetaMax,100);
        xGrid = cosd(angleGrid);
        yGrid = sind(angleGrid);
        circleMesh = linspace(rMin,rMax,Ncircles);
        spokeMesh = linspace(thetaMin,thetaMax,Nspokes);
        contour = abs((circleMesh - circleMesh(1))/Rrange+R(1)/Rrange);
        % plot circles
        for kk=1:length(contour)
            plot(xGrid*contour(kk), yGrid*contour(kk),'k:');
        end
        % radius tick label
        for kk=1:Ncircles
            
            position = 0.51.*(spokeMesh(min(Nspokes,round(Ncircles/2)))+...
                spokeMesh(min(Nspokes,1+round(Ncircles/2))));
            
            if abs(round(position)) ==90,
                % radial graduations
                text((contour(kk)).*cosd(90-position),...
                    (0.1+contour(kk)).*sind(86-position),...
                    num2str(circleMesh(kk),2),'verticalalignment','BaseLine',...
                    'horizontalAlignment', 'center',...
                    'handlevisibility','off','parent',cax);
                
                % annotate spokes
                text(contour(end).*0.6.*cosd(90-position),...
                    0.07+contour(end).*0.6.*sind(90-position),...
                    [labelR],'verticalalignment','bottom',...
                    'horizontalAlignment', 'right',...
                    'handlevisibility','off','parent',cax);
            else
                % radial graduations
                text((contour(kk)).*cosd(90-position),...
                    (contour(kk)).*sind(90-position),...
                    num2str(circleMesh(kk),2),'verticalalignment','BaseLine',...
                    'horizontalAlignment', 'right',...
                    'handlevisibility','off','parent',cax);
                
                % annotate spokes
                text(contour(end).*0.6.*cosd(90-position),...
                    contour(end).*0.6.*sind(90-position),...
                    [labelR],'verticalalignment','bottom',...
                    'horizontalAlignment', 'right',...
                    'handlevisibility','off','parent',cax);
            end
        end
        
    end
end

