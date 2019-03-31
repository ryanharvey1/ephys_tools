function analyze_cfc(Group1,Group2,VarName)


shadded_line_fig=figure; shadded_line_fig.Color=[1 1 1];shadded_line_fig.OuterPosition=[680 630 729 421];

% CONTROL

y=[nanmean(Group1(:,1)),nanmean(Group1(:,2)),nanmean(Group1(:,3)),nanmean(Group1(:,4))];

SEM=[nanstd(Group1(:,1))/sqrt(size(Group1(:,1),1)),nanstd(Group1(:,2))/sqrt(size(Group1(:,2),1)),...
    nanstd(Group1(:,3))/sqrt(size(Group1(:,3),1)),nanstd(Group1(:,4))/sqrt(size(Group1(:,4),1))];

h1=shadedErrorBar(1:4,y,SEM,'-k',0);
ylabel([VarName,'Power (dB/Hz)']);
ax1=gca; set(ax1,'XTick',[1 2 3 4],'Box','off','FontSize',20,'FontWeight','bold','LineWidth',2)

hold on

% PAE
y=[nanmean(Group2(:,1)),nanmean(Group2(:,2)),nanmean(Group2(:,3)),nanmean(Group2(:,4))];

SEM=[nanstd(Group2(:,1))/sqrt(size(Group2(:,1),1)),nanstd(Group2(:,2))/sqrt(size(Group2(:,2),1)),...
    nanstd(Group2(:,3))/sqrt(size(Group2(:,3),1)),nanstd(Group2(:,4))/sqrt(size(Group2(:,4),1))];

h2=shadedErrorBar(1:4,y,SEM,'-r',1);
set(h2.mainLine,'Color',[1 0 0])

% STATS

groups=[repmat({'Group1'},size(Group1,1),1);repmat({'Group2'},size(Group2,1),1)];
%     groups=data.textdata.Session1; groups(1,:)=[]; groups=groups(~cellfun(@isempty,groups));
measureVal=[Group1(:,1),Group1(:,2),Group1(:,3),Group1(:,4);Group2(:,1),Group2(:,2),Group2(:,3),Group2(:,4)];
t=table(groups,measureVal(:,1),measureVal(:,2),measureVal(:,3),measureVal(:,4),'VariableNames',{'species','meas1','meas2','meas3','meas4'});
Meas = table([1 2 3 4]','VariableNames',{'Measurements'});

% CREATE MODELS
% OVER ALL SESSIONS
rm1 = fitrm(t,'meas1-meas4~species','WithinDesign',Meas);

disp(strcat('--------------------------------------------------Quadrants--------------------------------------------------'))
disp('Sphericity')
mauchly(rm1)
disp(['Homogenity of Variance, p: ',num2str(vartestn(measureVal,'TestType','LeveneQuadratic','Display','off'))])
disp('Multivariate Tests *ALL SESSIONS*')
manovatbl=manova(rm1)
disp('Tests of Between-Subjects Effects *ALL SESSIONS*')
anovatbl=anova(rm1)


% MEAN AND SEM
disp(['Control  S1: ',num2str(nanmean(Group1(:,1))),' ± ',num2str(nanstd(Group1(:,1)/sqrt(length(Group1(:,1))))),...
    ', S2: ',num2str(nanmean(Group1(:,2))),' ± ',num2str(nanstd(Group1(:,2)/sqrt(length(Group1(:,2))))),...
    ', S3: ',num2str(nanmean(Group1(:,3))),' ± ',num2str(nanstd(Group1(:,3)/sqrt(length(Group1(:,3))))),...
    ', S4: ',num2str(nanmean(Group1(:,4))),' ± ',num2str(nanstd(Group1(:,4)/sqrt(length(Group1(:,4)))))])

disp(['PAE  S1: ',num2str(nanmean(Group2(:,1))),' ± ',num2str(nanstd(Group2(:,1)/sqrt(length(Group2(:,1))))),...
    ', S2: ',num2str(nanmean(Group2(:,2))),' ± ',num2str(nanstd(Group2(:,2)/sqrt(length(Group2(:,2))))),...
    ', S3: ',num2str(nanmean(Group2(:,3))),' ± ',num2str(nanstd(Group2(:,3)/sqrt(length(Group2(:,3))))),...
    ', S4: ',num2str(nanmean(Group2(:,4))),' ± ',num2str(nanstd(Group2(:,4)/sqrt(length(Group2(:,4)))))])

end