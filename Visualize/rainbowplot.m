function h=rainbowplot(y,y_minimum)
% rainbowplot: creates an area plot that is color coded by magnitude  
% 
% Input:
%           y           single vector of data
%           y_minimum   the minimum y value (default: 0)
% Output:
%           h           figure handle
% 
% Example:
%           y=rand(1,100);
%           rainbowplot(y,0);
% 
% Ryan E Harvey 2018

if isrow(y)
    y=y';
end
if ~exist('y_minimum','var')
    y_minimum=0;
end

h=fill([[1:length(y)]';length(y);1],[y;y_minimum;y_minimum],[y;y_minimum;y_minimum]);
xlim([1 length(y)])
ylim([y_minimum inf])
shading interp
% colormap jet
box off
end

