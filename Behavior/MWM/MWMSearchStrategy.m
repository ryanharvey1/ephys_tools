function Data=MWMSearchStrategy(path,radius,platformcenter)
% MWM SEARCH STRATEGY ANALYSIS
% 
% WISHLIST: 
%   CURVATURE (JARLIER ET AL. 2013)
% 
%   (RUEDIGER,2012) also defined in Illouz (2016)
%   THIGMOTAXIS 
%   SCANNING
%   DIRECTED SEARCH
%   RANDOM SWIM
%   CHAINING
%   FOCAL SEARCH
%   DIRECT SWIM

% Gehring search strategies
% TT,thigmotaxis
% IC,incursion
% SC,scanning,
% FS,focused search
% CR,chaining response
% SO,self orienting
% SS,scanning surroundings
% ST,target scanning
% DF,direct finding
% AT,approaching target

%   TURNING BEHAVIOR (HARVEY, 2009)
% 
% Ryan Harvey 2017
% 
% COMPUTE PERIMETER
ang=0:0.01:2*pi; perimeter=[radius*cos(ang)',radius*sin(ang)'];

% COMPUTE PLATFORM CORNERS 
TL=[platformcenter(1)-8,platformcenter(2)+8]; TR=[platformcenter(1)+8,platformcenter(2)+8];
BL=[platformcenter(1)-8,platformcenter(2)-8]; BR=[platformcenter(1)+8,platformcenter(2)-8];

% DIRECT 
x=linspace(path(1,1),platformcenter(1),100);
y=linspace(path(1,2),platformcenter(2),100);
right=x+8;
left=x-8;

% SPLIT UP PATH INTO SEGMENTS WITH 70% OVERLAP
% start to middle of platform
S2P=sqrt((platformcenter(2)-path(1,1))^2+(platformcenter(2)-path(1,2))^2);
% Find cumsum to split path into segments
CumSum=cumsum(sqrt(diff(path(:,1)).^2+diff(path(:,2)).^2));

for i=1:length(path)
    if i==1; SegIndex=S2P; else SegIndex=S2P+S2P*.3; end
    segment=path(CumSum<SegIndex,:);
    
    % DIRECT PATH
    
end

% SEARCH STRATEGY DECISION TREE





end