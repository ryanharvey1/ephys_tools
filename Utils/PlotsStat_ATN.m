function done=PlotsStat(Group1,Group2,VarNames)
for c=1:length(VarNames)
    shadded_line_fig=figure; shadded_line_fig.Color=[1 1 1];shadded_line_fig.OuterPosition=[680 630 729 421];
    % CONTROL
    y=[nanmean(Group1(:,c,1)),nanmean(Group1(:,c,2)),nanmean(Group1(:,c,3)),nanmean(Group1(:,c,4))];
    SEM=[nanstd(Group1(:,c,1))/sqrt(size(Group1(:,c,1),1)),nanstd(Group1(:,c,2))/sqrt(size(Group1(:,c,2),1)),...
        nanstd(Group1(:,c,3))/sqrt(size(Group1(:,c,3),1)),nanstd(Group1(:,c,4))/sqrt(size(Group1(:,c,4),1))];
    h1=shadedErrorBar(1:4,y,SEM,'-k',0);
    ylabel(VarNames(c));
    ax1=gca; set(ax1,'XTick',[1 2 3 4],'Box','off','FontSize',20,'FontWeight','bold','LineWidth',2)
    
    hold on
    
    % TILTED
    y=[nanmean(Group2(:,c,1)),nanmean(Group2(:,c,2)),nanmean(Group2(:,c,3)),nanmean(Group2(:,c,4))];
    SEM=[nanstd(Group2(:,c,1))/sqrt(size(Group2(:,c,1),1)),nanstd(Group2(:,c,2))/sqrt(size(Group2(:,c,2),1)),...
        nanstd(Group2(:,c,3))/sqrt(size(Group2(:,c,3),1)),nanstd(Group2(:,c,4))/sqrt(size(Group2(:,c,4),1))];
    h2=shadedErrorBar(1:4,y,SEM,'-r',1);
    set(h2.mainLine,'Color',[1 0 0])
    
    %     print(shadded_line_fig, '-dpdf', '-r300',{strcat(FigureLocation,filesep,VarNames(c),'_shadded_line_fig')})
    %     print(strcat(FigureLocation,filesep,VarNames(c),'_shadded_line_fig.pdf'),'-dpdf','-r300','-bestfit')
    %     print(shadded_line_fig,'-bestfit', '-dpdf', '-r600',char(strcat(FigureLocation,filesep,VarNames(c),'_shadded_line_fig')))
    
    %     close all
end

% STATS
for i=1:length(VarNames)
    groups=[repmat({'Group1'},size(Group1,1),1);repmat({'Group2'},size(Group2,1),1)];
%     groups=data.textdata.Session1; groups(1,:)=[]; groups=groups(~cellfun(@isempty,groups));
    measureVal=[Group1(:,i,1),Group1(:,i,2),Group1(:,i,3),Group1(:,i,4);Group2(:,i,1),Group2(:,i,2),Group2(:,i,3),Group2(:,i,4)];
    t=table(groups,measureVal(:,1),measureVal(:,2),measureVal(:,3),measureVal(:,4),'VariableNames',{'species','meas1','meas2','meas3','meas4'});
    Meas = table([1 2 3 4]','VariableNames',{'Measurements'});
    
    % CREATE MODELS
    % OVER ALL SESSIONS
    rm1 = fitrm(t,'meas1-meas4~species','WithinDesign',Meas);
    % SESSIONS 1 2 3 - rotation
    rm2 = fitrm(t,'meas2-meas4~species','WithinDesign',Meas(2:4,1));
    % SESSIONS 1 4 5 - dark
    rm3 = fitrm(t,'meas1,meas2,meas4~species','WithinDesign',Meas([1,2,4],1));
    % SESSIONS 1 3 5 - stability of standard sessions
    rm4 = fitrm(t,'meas1,meas2~species','WithinDesign',Meas([1,2],1));
    
    % tbl = mauchly(rm)
    % ranovatbl = ranova(rm)
    disp(strcat('--------------------------------------------------',VarNames(i),'--------------------------------------------------'))
    disp('Sphericity')
    mauchly(rm1)
    disp(['Homogenity of Variance, p: ',num2str(vartestn(measureVal,'TestType','LeveneQuadratic','Display','off'))])
    disp('Multivariate Tests *ALL SESSIONS*')
    manovatbl=manova(rm1)
    disp('Tests of Between-Subjects Effects *ALL SESSIONS*')
    anovatbl=anova(rm1)
    
    disp('Multivariate Tests *ROTATION*')
    disp('Sphericity')
    mauchly(rm2)
    disp(['Homogenity of Variance, p: ',num2str(vartestn(measureVal(:,2:4),'TestType','LeveneQuadratic','Display','off'))])
    manovatbl=manova(rm2)
    disp('Tests of Between-Subjects Effects *ROTATION*')
    anova(rm2)
    
    disp('Multivariate Tests *DARK*')
    disp('Sphericity')
    mauchly(rm3)
    disp(['Homogenity of Variance, p: ',num2str(vartestn(measureVal(:,[1,2,4]),'TestType','LeveneQuadratic','Display','off'))])
    manovatbl=manova(rm3)
    disp('Tests of Between-Subjects Effects *DARK*')
    anova(rm3)
    
    disp('Multivariate Tests *STABILITY*')
    disp('Sphericity')
    mauchly(rm4)
    disp(['Homogenity of Variance, p: ',num2str(vartestn(measureVal(:,[1,2]),'TestType','LeveneQuadratic','Display','off'))])
    manovatbl=manova(rm4)
    disp('Tests of Between-Subjects Effects *STABILITY*')
    anova(rm4)
end

% MEAN AND SEM
for v=1:length(VarNames)
    disp(VarNames(v))
    disp(['Control  S1: ',num2str(nanmean(Group1(:,v,1))),' ± ',num2str(nanstd(Group1(:,v,1)/sqrt(length(Group1(:,v,1))))),...
        ', S2: ',num2str(nanmean(Group1(:,v,2))),' ± ',num2str(nanstd(Group1(:,v,2)/sqrt(length(Group1(:,v,2))))),...
        ', S3: ',num2str(nanmean(Group1(:,v,3))),' ± ',num2str(nanstd(Group1(:,v,3)/sqrt(length(Group1(:,v,3))))),...
        ', S4: ',num2str(nanmean(Group1(:,v,4))),' ± ',num2str(nanstd(Group1(:,v,4)/sqrt(length(Group1(:,v,4))))),...
        ', S5: ',num2str(nanmean(Group1(:,v,5))),' ± ',num2str(nanstd(Group1(:,v,5)/sqrt(length(Group1(:,v,5)))))])
    
    disp(['PAE  S1: ',num2str(nanmean(Group2(:,v,1))),' ± ',num2str(nanstd(Group2(:,v,1)/sqrt(length(Group2(:,v,1))))),...
        ', S2: ',num2str(nanmean(Group2(:,v,2))),' ± ',num2str(nanstd(Group2(:,v,2)/sqrt(length(Group2(:,v,2))))),...
        ', S3: ',num2str(nanmean(Group2(:,v,3))),' ± ',num2str(nanstd(Group2(:,v,3)/sqrt(length(Group2(:,v,3))))),...
        ', S4: ',num2str(nanmean(Group2(:,v,4))),' ± ',num2str(nanstd(Group2(:,v,4)/sqrt(length(Group2(:,v,4))))),...
        ', S5: ',num2str(nanmean(Group2(:,v,5))),' ± ',num2str(nanstd(Group2(:,v,5)/sqrt(length(Group2(:,v,5)))))])
end
done=1;
end