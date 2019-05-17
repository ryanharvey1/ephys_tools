function plotspread_wrapper(x1,x2)
% plotspread_wrapper: wrapper for plotSpread
% Takes in data from 2 groups and creates bee swarm plots
%
%   Input: 
%           x1,x2: two vectors of data to plot
%
% Ryan Harvey 2019

data = {x1, x2};

% check for plotSpread package
if exist('plotSpread','file')~=2
    error('Add plotSpread to path')
end

plotSpread(data, ...
    'xNames', {'control','lesion'}, ...
    'distributionMarkers', {'.', '.'},...
    'distributionColors',{'k','r'});

hold on
plot([.5,1.5],[nanmean(data{1}),nanmean(data{1})],'Color','k','LineWidth',2)
plot([1.5,2.5],[nanmean(data{2}),nanmean(data{2})],'Color','r','LineWidth',2)

sem=nanstd(data{1})/sqrt(length(data{1}));
v = [.5 nanmean(data{1})-sem; 1.5 nanmean(data{1})-sem; 1.5 nanmean(data{1})+sem; 0.5 nanmean(data{1})+sem];
f = [1 2 3 4];
patch('Faces',f,'Vertices',v,'FaceColor','black','FaceAlpha',.3,'EdgeColor','none')

sem=nanstd(data{2})/sqrt(length(data{2}));
v = [1.5 nanmean(data{2})-sem; 2.5 nanmean(data{2})-sem; 2.5 nanmean(data{2})+sem; 1.5 nanmean(data{2})+sem];
f = [1 2 3 4];
patch('Faces',f,'Vertices',v,'FaceColor','red','FaceAlpha',.3,'EdgeColor','none')

chi=get(gca, 'Children');
set(gca, 'Children',flipud(chi),'LineWidth',1,'FontSize',12,'FontWeight','bold');
end
