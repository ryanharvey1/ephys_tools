% analyzeigoroutput
clear;clc;close all
addpath('/Users/RyanHarvey/GoogleDrive/MatlabDir/BClarkToolbox/Analysis')
load('/Users/ryanharvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/Ryan_et_al_2018_MasterData.mat')


%%
% load('/Users/RyanHarvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/NewPipelineResults.mat')
% % [PeakRate,nSpikes,OverallFR,NumbActiveBins,sparsity,InformationContent,Coherence,Field2Wall,borderScore,FieldWidth,infieldFR,outfieldFR,E,c,p];
%
% results.control.results(i,1,:)
%
% [PeakRate,nSpikes,OverallFR,NumbActiveBins,sparsity,InformationContent,Coherence,Field2Wall,borderScore,FieldWidth,infieldFR,outfieldFR,E,c,p];

% control
% fig 1 example cells
% %  find firing rate values for place cell examples
% controlexamples=[3,5,6,11,12,17,34];
% for s=1:5
%     for ii=1:length(controlexamples)
%         controlmeasures(ii,s)=tempdatcontrol(controlexamples(ii),s);
%     end
% end
%
%
% tiltedexamples=[19,38,43,68,69,75,77];
% for s=1:5
%     for ii=1:length(tiltedexamples)
%         tiltedmeasures(ii,s)=tempdattilted(tiltedexamples(ii),s);
%     end
% end

% for i=1:length(tiltedexamples)
%     tiltedmeasures{i}=results.tilted.results(tiltedexamples(i),1,:);
% end
%% PLOT RAW DATA
controlexamples=[3,5,6,11,12,17,34];

tiltedexamples=[19,38,43,68,69,75,77];

for i=[6,13,2]
    tempdatcontrol=[controlplacecelldata(:,i,1),controlplacecelldata(:,i,2),controlplacecelldata(:,i,3),controlplacecelldata(:,i,4),controlplacecelldata(:,i,5)];
    tempdattilted=[tiltedplacecelldata(:,i,1),tiltedplacecelldata(:,i,2),tiltedplacecelldata(:,i,3),tiltedplacecelldata(:,i,4),tiltedplacecelldata(:,i,5)];
    for s=1:5
        cdf_fig=figure; cdf_fig.Color=[1 1 1]; cdf_fig.Position=[768 259 401 431];
        [f1,x1] = ecdf(tempdatcontrol(:,s));
        [f2,x2] = ecdf(tempdattilted(:,s));
        
        p1=plot(x1,f1);
        set(p1,'LineWidth',4,'Color','k')
        hold on
        p2=plot(x2,f2);
        set(p2,'LineWidth',4,'Color','r')
        ylabel('Cumulative Frequency')
        xlabel(Varnames(i))
        title(['Session',num2str(s)])
        set(gca,'FontSize',20,'FontWeight','bold','LineWidth',2,'box','off')
        
        % SCATTER EXAMPLE CELLS
        % control
        for ii=1:length(controlexamples)
            controlmeasures(ii,1)=tempdatcontrol(controlexamples(ii),s);
        end
        [C,ia,ib]=intersect(x1,controlmeasures);
        h1=scatter(x1(ia),f1(ia),'filled','MarkerFaceColor',[.75 .75 .75],'MarkerEdgeColor','k');
        set(h1,'SizeData',250);
        
        % tilted
        for ii=1:length(tiltedexamples)
            tiltedmeasures(ii,1)=tempdattilted(tiltedexamples(ii),s);
        end
        [C,ia,ib]=intersect(x2,tiltedmeasures);
        h2=scatter(x2(ia),f2(ia),'filled','MarkerFaceColor',[.75 .75 .75],'MarkerEdgeColor','k');
        set(h2,'SizeData',250);
        
        
        print(cdf_fig,'-bestfit', '-dpdf', '-r600',['/Users/RyanHarvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/Figures_For_Paper/Example Cells on CDF',filesep,'cdf_Fig','_',Varnames{i},'_','Session',num2str(s),'.pdf'])
        close all
    end
    
end


% handles = plotSpread(tempdattilted,'distributionMarkers','.','distributionColors',{'r','r','r','r','r'})
%
% handles = plotSpread(tempdatcontrol,'categoryLabels',{'1','2','3','4','5'},'distributionMarkers','o','distributionColors',{'k','k','k','k','k'});hold on

% figure;
% y=sort(tempdatcontrol(:,1),'descend');
% scatter(linspace(.5,1.5,length(tempdatcontrol)),y);hold on
%
% y=sort(tempdatcontrol(:,2),'descend');
% scatter(linspace(1.5,2.5,length(tempdatcontrol)),y);hold on
%
% y=sort(tempdatcontrol(:,3),'descend');
% scatter(linspace(2.5,3.5,length(tempdatcontrol)),y);hold on
%
% y=sort(tempdatcontrol(:,4),'descend');
% scatter(linspace(3.5,4.5,length(tempdatcontrol)),y);hold on
%
% y=sort(tempdatcontrol(:,5),'descend');
% scatter(linspace(4.5,5.5,length(tempdatcontrol)),y);hold on

%% SCATTER DISTRIBUTION FOR ALL CELLS AND PLACE CELLS

% plot all cells
scatter_fig=figure; scatter_fig.Color=[1 1 1];
plot(controldata(:,3,1),controldata(:,2,1),'.k');hold on
plot(tilteddata(:,3,1),tilteddata(:,2,1),'.r');
xlabel('Information Content (bits/spk)')
ylabel('Coherence')
set(gca,'Box','off','FontWeight','bold','FontSize',15,'LineWidth',3)

% plot cut off lines
l1=plot([.6;.6],[0;1],'k');
l2=plot([0;max([tilteddata(:,3,1);controldata(:,3,1)])],[.5;.5],'k');
set(l1,'LineWidth',3)
set(l2,'LineWidth',3)

% plot place cells
scatter(controlplacecelldata(:,3,1),controlplacecelldata(:,2,1),40,'k','filled');hold on
scatter(tiltedplacecelldata(:,3,1),tiltedplacecelldata(:,2,1),40,'r','filled');

controlplacecellpaths(controlplacecelldata(:,3,1)<.5)
controlplacecelldata(controlplacecelldata(:,3,1)<.5,3)

tiltedplacecellpaths(tiltedplacecelldata(:,2,1)<.5)
tiltedplacecelldata(tiltedplacecelldata(:,2,1)<.5,2)


%% SCATTER DISTRIBUTION FOR PRINCIPLE CELLS AND PLACE CELLS

% filter down to just principle cells
controlfilt=controldata(controldata(:,4,1)<10,:,1);
controlfilt=controlfilt(controlfilt(:,5)>185 | controlfilt(:,5)==0 | isnan(controlfilt(:,5)),:);

tiltedfilt=tilteddata(tilteddata(:,4,1)<10,:,1);
tiltedfilt=tiltedfilt(tiltedfilt(:,5)>185 | tiltedfilt(:,5)==0 | isnan(tiltedfilt(:,5)),:);

% plot cells
scatter_fig=figure; scatter_fig.Color=[1 1 1];
scatter(controlfilt(:,3,1),controlfilt(:,2,1),'k','filled');hold on
scatter(tiltedfilt(:,3,1),tiltedfilt(:,2,1),'r','filled');
xlabel('Information Content (bits/spk)')
ylabel('Coherence')
set(gca,'Box','off','FontWeight','bold','FontSize',15,'LineWidth',3)

% % plot cut off lines
% l1=plot([.6;.6],[0;1],'k');
% l2=plot([0;max([tilteddata(:,3,1);controldata(:,3,1)])],[.5;.5],'k');
% set(l1,'LineWidth',3)
% set(l2,'LineWidth',3)

% % plot place cells
% scatter(controlplacecelldata(:,3,1),controlplacecelldata(:,2,1),40,'k','filled');hold on
% scatter(tiltedplacecelldata(:,3,1),tiltedplacecelldata(:,2,1),40,'r','filled');


%% cdf of all principle cells and each place cells position

% filter down to just principle cells
controlfilt=controldata(controldata(:,4,1)<10 & controldata(:,5)>185 | controldata(:,5)==0 | isnan(controldata(:,5)),:,1);
controlfiltpaths=allcontrolpaths(controldata(:,4,1)<10 & controldata(:,5)>185 | controldata(:,5)==0 | isnan(controldata(:,5)),1,1);

tiltedfilt=tilteddata(tilteddata(:,4,1)<10 & tilteddata(:,5)>185 | tilteddata(:,5)==0 | isnan(tilteddata(:,5)),:,1);
tiltedfiltpaths=alltiltedpaths(tilteddata(:,4,1)<10 & tilteddata(:,5)>185 | tilteddata(:,5)==0 | isnan(tilteddata(:,5)),1,1);



% [ AllStats ] = CDFplots(controlfilt(:,2:3),tiltedfilt(:,2:3),{'Control','Tilted'},{'Coherence','Infocontent'},1)

% filter out place cells
controlNoPlaceCell=controlfilt(~contains(controlfiltpaths(:,1,1),controlplacecellpaths),:,1);
tiltedNoPlaceCell=tiltedfilt(~contains(tiltedfiltpaths(:,1,1),tiltedplacecellpaths),:,1);


for i=2:3
    [ AllStats ] = CDFplots(controlNoPlaceCell(:,i),tiltedNoPlaceCell(:,i),{'Control','Tilted'},Varnames{i},2)
        cdf_fig=figure(1); cdf_fig.Color=[1 1 1]; cdf_fig.Position=[768 259 401 431]; 

        print(cdf_fig,'-bestfit', '-dpdf', '-r600',...
        ['/Users/RyanHarvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/Figures_For_Paper/PlaceCellsOverAllPrincipleCells',...
        filesep,'noplacecells','_',Varnames{i},'.pdf'])
    close all
end



for i=1:size(controlfilt,2)
    cdf_fig=figure; cdf_fig.Color=[1 1 1]; cdf_fig.Position=[768 259 401 431];
    [f1,x1] = ecdf(controlfilt(:,i));
    [f2,x2] = ecdf(tiltedfilt(:,i));
    
    p1=plot(x1,f1);
    set(p1,'LineWidth',4,'Color','k')
    hold on
    p2=plot(x2,f2);
    set(p2,'LineWidth',4,'Color','r')
    ylabel('Cumulative Frequency')
    xlabel(Varnames(i))
    set(gca,'FontSize',20,'FontWeight','bold','LineWidth',2,'box','off')
    
%     % SCATTER EXAMPLE CELLS
%     % control
% %     for ii=1:length(controlexamples)
% %         controlmeasures(ii,1)=tempdatcontrol(controlexamples(ii),i);
% %     end
%     [C,ia,ib]=intersect(x1,controlplacecelldata(:,i,1));
%     h1=scatter(x1(ia),f1(ia),'filled','MarkerFaceColor',[0 0 0],'MarkerEdgeColor','k');
%     set(h1,'SizeData',100);
%     
%     % tilted
% %     for ii=1:length(tiltedexamples)
% %         tiltedmeasures(ii,1)=tempdattilted(tiltedexamples(ii),i);
% %     end
%     [C,ia,ib]=intersect(x2,tiltedplacecelldata(:,i,1));
%     h2=scatter(x2(ia),f2(ia),'filled','MarkerFaceColor',[1 0 0],'MarkerEdgeColor','k');
%     set(h2,'SizeData',100);
    
    
    print(cdf_fig,'-bestfit', '-dpdf', '-r600',...
        ['/Users/RyanHarvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/Figures_For_Paper/PlaceCellsOverAllPrincipleCells',...
        filesep,'cdf_Fig','_',Varnames{i},'.pdf'])
    close all
end

%% histogram of all principle cells and each place cells position


% filter down to just principle cells
controlfilt=controldata(controldata(:,4,1)<10,:,1);
controlfilt=controlfilt(controlfilt(:,5)>185 | controlfilt(:,5)==0 | isnan(controlfilt(:,5)),:);

tiltedfilt=tilteddata(tilteddata(:,4,1)<10,:,1);
tiltedfilt=tiltedfilt(tiltedfilt(:,5)>185 | tiltedfilt(:,5)==0 | isnan(tiltedfilt(:,5)),:);

% find max min
edges=[linspace(min([controlfilt(:,i);tiltedfilt(:,i)]),max([controlfilt(:,i);tiltedfilt(:,i)]),100)];

% bin
controlbinded=histcounts(controlfilt,edges);
tiltedbinded=histcounts(tiltedfilt,edges);

% rescale 0 to 1
controlbinded=rescale(controlbinded,0,1);
tiltedbinded=rescale(tiltedbinded,0,1);

% plot
figure;
plot(edges(2:end),controlbinded,'k');hold on
plot(edges(2:end),tiltedbinded,'r');



controlbinded=histcounts(controlplacecelldata,edges);
tiltedbinded=histcounts(tiltedplacecelldata,edges);

controlbinded=rescale(controlbinded,0,1);
tiltedbinded=rescale(tiltedbinded,0,1);
figure;
area(edges(2:end),controlbinded);hold on
area(edges(2:end),tiltedbinded);




%% BAR PLOTS
Group1=controlplacecelldata;
Group2=tiltedplacecelldata;
for c=1:length(Varnames)
    shadded_line_fig=figure; shadded_line_fig.Color=[1 1 1];shadded_line_fig.OuterPosition=[680 630 729 421];
    % CONTROL
    y=[nanmean(Group1(:,c,1)),nanmean(Group1(:,c,2)),nanmean(Group1(:,c,3)),nanmean(Group1(:,c,4)),nanmean(Group1(:,c,5))];
    SEM=[nanstd(Group1(:,c,1))/sqrt(size(Group1(:,c,1),1)),nanstd(Group1(:,c,2))/sqrt(size(Group1(:,c,2),1)),...
        nanstd(Group1(:,c,3))/sqrt(size(Group1(:,c,3),1)),nanstd(Group1(:,c,4))/sqrt(size(Group1(:,c,4),1)),...
        nanstd(Group1(:,c,5))/sqrt(size(Group1(:,c,5),1))];
    h1=shadedErrorBar(1:5,y,SEM,'-k',0);
    ylabel(Varnames(c));
    ax1=gca; set(ax1,'XTick',[1 2 3 4 5],'Box','off','FontSize',20,'FontWeight','bold','LineWidth',2)
    
    hold on
    
    % TILTED
    y=[nanmean(Group2(:,c,1)),nanmean(Group2(:,c,2)),nanmean(Group2(:,c,3)),nanmean(Group2(:,c,4)),nanmean(Group2(:,c,5))];
    SEM=[nanstd(Group2(:,c,1))/sqrt(size(Group2(:,c,1),1)),nanstd(Group2(:,c,2))/sqrt(size(Group2(:,c,2),1)),...
        nanstd(Group2(:,c,3))/sqrt(size(Group2(:,c,3),1)),nanstd(Group2(:,c,4))/sqrt(size(Group2(:,c,4),1)),...
        nanstd(Group2(:,c,5))/sqrt(size(Group2(:,c,5),1))];
    h2=shadedErrorBar(1:5,y,SEM,'-r',1);
    set(h2.mainLine,'Color',[1 0 0])
    
    %     print(shadded_line_fig, '-dpdf', '-r300',{strcat(FigureLocation,filesep,VarNames(c),'_shadded_line_fig')})
    %     print(strcat(FigureLocation,filesep,VarNames(c),'_shadded_line_fig.pdf'),'-dpdf','-r300','-bestfit')
    %     print(shadded_line_fig,'-bestfit', '-dpdf', '-r600',char(strcat(FigureLocation,filesep,VarNames(c),'_shadded_line_fig')))
    
    %     close all
end
clear y shadded_line_fig SEM h1 h2 Group1 Group2 C ax1




% load('/Users/RyanHarvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/RHPlaceCellDataWorksheet_Compiled.mat')
%
%
%
% % load other non-igor data
% varnames=data.textdata.Field_StatsPlaceCells0x2DTilted(1,3:end);
% % extract control
% control=data.data.Field_StatsPlaceCells0x2DTilted(contains(data.textdata.Field_StatsPlaceCells0x2DTilted(2:end,2),'Control'),:);
% % extract tilted
% tilted=data.data.Field_StatsPlaceCells0x2DTilted(contains(data.textdata.Field_StatsPlaceCells0x2DTilted(2:end,2),'Tilted'),:);
%
%
% controlpaths=data.textdata.Field_StatsPlaceCells0x2DTilted(contains(data.textdata.Field_StatsPlaceCells0x2DTilted(:,2),'Control'),:);
% tiltedpaths=data.textdata.Field_StatsPlaceCells0x2DTilted(contains(data.textdata.Field_StatsPlaceCells0x2DTilted(:,2),'Tilted'),:);
%
%
% load('/Users/RyanHarvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/NewData.mat')
%
% % contains(data.textdata.Field_StatsPlaceCells0x2DTilted(2:end,1),newdata.textdata(2:end,1));
% placecells=data.textdata.Field_StatsPlaceCells0x2DTilted(2:end,1);
%
% varnames2=newdata.textdata(1,14:17);
%
% control=[control,newdata.data(contains(newdata.textdata(2:end,1),controlpaths(:,1),'IgnoreCase',true),13:16)];
% tilted=[tilted,newdata.data(contains(newdata.textdata(2:end,1),tiltedpaths(:,1),'IgnoreCase',true),13:16)];
%
%
% load('/Users/RyanHarvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/NewData2.mat');
% newdata2=ans;
%
% varnames=[varnames,varnames2,{'infieldFR','outfieldFR'}];
%
% control=[control,newdata2.data(contains(newdata.textdata(2:end,1),controlpaths(:,1),'IgnoreCase',true),17:18)];
% tilted=[tilted,newdata2.data(contains(newdata.textdata(2:end,1),tiltedpaths(:,1),'IgnoreCase',true),17:18)];
%
%
% controlmain=control;
% tiltedmain=tilted;
%
% control=control(:,[1:3,7,10,12:end]);
% tilted=tilted(:,[1:3,7,10,12:end]);
% varnamesreduced=varnames(:,[1:3,7,10,12:end]);
%
% control3d(:,:,1)=[controlplacecelldata(:,:,1),control(1:5:end,:)];
% control3d(:,:,2)=[controlplacecelldata(:,:,2),control(2:5:end,:)];
% control3d(:,:,3)=[controlplacecelldata(:,:,3),control(3:5:end,:)];
% control3d(:,:,4)=[controlplacecelldata(:,:,4),control(4:5:end,:)];
% control3d(:,:,5)=[controlplacecelldata(:,:,5),control(5:5:end,:)];
%
% tilted3d(:,:,1)=[tiltedplacecelldata(:,:,1),tilted(1:5:end,:)];
% tilted3d(:,:,2)=[tiltedplacecelldata(:,:,2),tilted(2:5:end,:)];
% tilted3d(:,:,3)=[tiltedplacecelldata(:,:,3),tilted(3:5:end,:)];
% tilted3d(:,:,4)=[tiltedplacecelldata(:,:,4),tilted(4:5:end,:)];
% tilted3d(:,:,5)=[tiltedplacecelldata(:,:,5),tilted(5:5:end,:)];
%
%
% Varnames=[{'PeakRateIGOR','avgrateIGOR','coherenceIGOR','infocontentIGOR','APdurationIGOR'},varnamesreduced];
%
% placecellresults.control(:,[1:3,7,10,12,13:18]);
% placecellresults.varnames(:,[1:3,7,10,12,13:18])
%
% controlplacecelldata=control3d;
% tiltedplacecelldata=tilted3d;