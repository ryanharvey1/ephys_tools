% Cued_Opposite_Ben
close all
S=importdata('/Users/RyanHarvey/Downloads/Expt1_Speed.xlsx');
rats=fieldnames(S.data);

% PLOT EVERY LINE
for i=1:length(rats)
    temp=S.data.(rats{i});
    cued=temp(:,1);
    opposite=temp(:,2);

    cued(isnan(cued))=[];
    opposite(isnan(opposite))=[];
    
    sizes(i,:)=[length(cued),length(opposite)];
    
    subplot(1,2,1)
    figure(1);h1=plot(linspace(0,1,length(cued)),cued);title('Cued');hold on
    subplot(1,2,2)
    figure(1);h2=plot(linspace(0,1,length(opposite)),opposite);title('Opposite');hold on    
end

% PLOT MEAN
up=max(max(sizes));
for i=1:length(rats)
    temp=S.data.(rats{i});
    cued=temp(:,1);
    opposite=temp(:,2);

    cued(isnan(cued))=[];
    opposite(isnan(opposite))=[];
    
    cuedALL(:,i)=imresize(cued,[up,1]);
    oppositeALL(:,i)=imresize(opposite,[up,1]);  
end
cuedmean=mean(cuedALL,2);
oppositemean=mean(oppositeALL,2);


shadded_line_fig=figure; shadded_line_fig.Color=[1 1 1];shadded_line_fig.OuterPosition=[680 630 729 421];
h1=shadedErrorBar(linspace(0,1,length(cuedmean)),cuedmean,std(cuedALL,[],2)/sqrt(size(cuedALL,2)),'-k',0);
ylabel('Velocity');
format shortg;
ax1=gca; set(ax1,'XTick',round(linspace(0,1,11),1,'significant'),'Box','off','FontSize',20,'FontWeight','bold','LineWidth',2)
hold on
h2=shadedErrorBar(linspace(0,1,length(oppositemean)),oppositemean,std(oppositeALL,[],2)/sqrt(size(oppositeALL,2)),'-r',1);
set(h2.mainLine,'Color',[1 0 0])

print(shadded_line_fig,'-bestfit', '-dpdf', '-r600',char(strcat('/Users/RyanHarvey/Downloads',filesep,'Velo','_shadded_line_fig')))

