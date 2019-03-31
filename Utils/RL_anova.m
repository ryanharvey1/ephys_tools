function RL_anova(Group1,Group2,VarName)

shadded_line_fig=figure; shadded_line_fig.Color=[1 1 1];shadded_line_fig.OuterPosition=[680 630 729 421];

% CONTROL
SEM=nanstd(Group1)/sqrt(size(Group1,1));
h1=shadedErrorBar(1:size(Group1,2),nanmean(Group1),SEM,'-k',0);
ylabel(VarName);
ax1=gca; set(ax1,'XTick',1:size(Group1,2),'Box','off','FontSize',20,'FontWeight','bold','LineWidth',2)
xlim([1 size(Group1,2)])

hold on

SEM=nanstd(Group2)/sqrt(size(Group2,1));
h2=shadedErrorBar(1:size(Group2,2),nanmean(Group2),SEM,'-r',1);
set(h2.mainLine,'Color',[1 0 0])

% STATS

groups=[repmat({'Group1'},size(Group1,1),1);repmat({'Group2'},size(Group2,1),1)];
%     groups=data.textdata.Session1; groups(1,:)=[]; groups=groups(~cellfun(@isempty,groups));
measureVal=[Group1;Group2];

repmat({'meas'},1,size(Group1,2))
for i=1:size(Group1,2)
    measures{i}=['meas',num2str(i)];  
end

t=array2table([groups,num2cell(measureVal)],'VariableNames',['species',measures]);
for i=1:length(measures)
    t.(measures{i})=[t.(measures{i}){:}]';
end
Meas = table([1:size(Group1,2)]','VariableNames',{'Measurements'});

% CREATE MODELS
% OVER ALL SESSIONS
rm1 = fitrm(t,['meas1-',measures{end},' ~ species'],'WithinDesign',Meas);


disp(strcat('--------------------------------------------------',VarName,'--------------------------------------------------'))
disp('Sphericity')
mauchly(rm1)
disp(['Homogenity of Variance, p: ',num2str(vartestn(measureVal,'TestType','LeveneQuadratic','Display','off'))])
disp('Multivariate Tests *ALL SESSIONS*')
manovatbl=manova(rm1)
disp('Tests of Between-Subjects Effects *ALL SESSIONS*')
anovatbl=anova(rm1)


% MEAN AND SEM
% disp(['Control  S1: ',num2str(nanmean(Group1(:,1))),' ± ',num2str(nanstd(Group1(:,1)/sqrt(length(Group1(:,1))))),...
%     ', S2: ',num2str(nanmean(Group1(:,2))),' ± ',num2str(nanstd(Group1(:,2)/sqrt(length(Group1(:,2))))),...
%     ', S3: ',num2str(nanmean(Group1(:,3))),' ± ',num2str(nanstd(Group1(:,3)/sqrt(length(Group1(:,3))))),...
%     ', S4: ',num2str(nanmean(Group1(:,4))),' ± ',num2str(nanstd(Group1(:,4)/sqrt(length(Group1(:,4)))))])
% 
% disp(['PAE  S1: ',num2str(nanmean(Group2(:,1))),' ± ',num2str(nanstd(Group2(:,1)/sqrt(length(Group2(:,1))))),...
%     ', S2: ',num2str(nanmean(Group2(:,2))),' ± ',num2str(nanstd(Group2(:,2)/sqrt(length(Group2(:,2))))),...
%     ', S3: ',num2str(nanmean(Group2(:,3))),' ± ',num2str(nanstd(Group2(:,3)/sqrt(length(Group2(:,3))))),...
%     ', S4: ',num2str(nanmean(Group2(:,4))),' ± ',num2str(nanstd(Group2(:,4)/sqrt(length(Group2(:,4)))))])

end