% MWMwithWMMML
clear; close all; clc
% load('/Users/RyanHarvey/Dropbox/school work/UNM/Lab/Projects/ClarkP30_WaterMaze/ProbeT2/Mclassification/class_242_1715_90_07_10_10_mr0-11Jul2017h15m26/merged_5.mat','-mat')
load('D:\ClarkP30 Tracking\Segmentation Analysis\T3\1RefMem\Analysis_91117\T3_RefMem_91117\Mclassification\class_403_4357_150_08_1_1_mr0-11Sep2017h16m25\merged_1.mat')
Labeled=classification_configs.CLASSIFICATION.class_map;
Features=classification_configs.FEATURES; 


% SPLIT UP BY SESSION
for i=1:size(classification_configs.CLASSIFICATION.segments.items,2)
    group=classification_configs.CLASSIFICATION.segments.items(1,i).group;
    ID=classification_configs.CLASSIFICATION.segments.items(1,i).id;
    session=classification_configs.CLASSIFICATION.segments.items(1,i).track;session(end-3:end)=[];
%     Strategy=Labeled(i);
%     Feature=Features(i,:);
    sessions(i,:)={session};
%     results(i,:)=[ID,group,Strategy,Feature];
    results(i,:)=[ID,group];
end
RESULTS=[sessions,num2cell(results),num2cell(Labeled'),num2cell(Features)];

% WRITE RESULTS
ResultTable=cell2table(RESULTS,'VariableNames',{'Session','AnimalID','Group','Strategy','MedianRadus',...
    'IQRRadius','Focus','CentreDisplacement','InnerRadiusVariation','PlatformProximity','BoundaryEccentricity','LongestLoop','Latency','PathLength','AverageSpeed'});
writetable(ResultTable,['D:\ClarkP30 Tracking\Segmentation Analysis\T3\1RefMem\Analysis_91117\T3_RefMem_91117' filesep '_results','.csv'],'WriteRowNames',true);


% PLOT ALL STRATEGIES
NumofStrat=unique(classification_configs.CLASSIFICATION.class_map);
for ii=0:length(NumofStrat)
    [I]=find(classification_configs.CLASSIFICATION.class_map==ii);
    for i=1:length(I)
        figure(ii+1)
        plot(classification_configs.CLASSIFICATION.segments.items(1, I(i)).points(:,2),classification_configs.CLASSIFICATION.segments.items(1, I(i)).points(:,3))
        hold on
    end
end
