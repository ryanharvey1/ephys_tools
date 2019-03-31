function [freq,p4]=ECDF_plot(data,varname)
% ECDF_plot plots the ecdf 
% 
% Input:        
%           Data:       structured array* of data from each group
%           varname:    string of variable name (ex. 'Latency')
% Output:
%           freq:       structured array of ecdf'd data for statistics 
% 
% ////////////////////////////////////////////////////////////////////////
% *Example of structured array
% 
% data.group1=group1(:,1);
% data.group2=group2(:,1);
% data.group3=group3(:,1);
% data.group4=group4(:,1);
% 
% \/ Further break down of data structure \/
% 
% data          .           group1          =           group1(:,1);
%  /\                         /\                            /\
% Main Structure           group name                  data from group
% ////////////////////////////////////////////////////////////////////////
% 
% Ryan E. Harvey 2017

groupname=fieldnames(data);
for i=1:size(groupname)
    [f,x] = ecdf(data.(groupname{i}));
    p4=plot(x,f); 
    freq.(groupname{i})=x;
    set(p4,'LineWidth',4)
    hold on
end
box off
legend(groupname,'FontSize',12,'Location','best','box','off','FontSize',15,'FontWeight','Bold')
xlabel(varname);ylabel('Cumulative Frequency')
ax=gca;
set(ax,'FontSize',20,'FontWeight','bold','LineWidth',2)
end