% AnalyzeResults_PAE_Cylinder
clear
data=compileResults('D:\Projects\PAE_PlaceCell\ProcessedData');

control={'RH13','RH14','LS21','LS23','LE2821','LE2823','LEM3116','LEM3120'};
pae={'RH11','RH16','LS17','LS19','LE2813','LEM3124'};

%% COMPILE GROUPS
data.control.measures=[];
data.control.id=[];
for i=1:length(control)
    data.control.measures=cat(1,data.control.measures,data.(control{i}).measures);
    data.control.id=cat(1,data.control.id,data.(control{i}).id);
end

data.pae.measures=[];
data.pae.id=[];
for i=1:length(pae)
    data.pae.measures=cat(1,data.pae.measures,data.(pae{i}).measures);
    data.pae.id=cat(1,data.pae.id,data.(pae{i}).id);
end

%% COMPILE Cylinder data
session=3;
group1=data.control.measures(:,:,session);
group2=data.pae.measures(:,:,session);
group1id=data.control.id;
group2id=data.pae.id;

%% Delete rows with all nans (sessions without cylinder)
to_delete=sum(isnan(group1),2)==size(group1,2);
group1(to_delete,:)=[];
group1id(to_delete,:)=[];

to_delete=sum(isnan(group2),2)==size(group2,2);
group2(to_delete,:)=[];
group2id(to_delete,:)=[];

%% DELETE MEASURES FOR OPEN ARENA
varnames=data.varnames;
% colstodelete=size(group1,1)==sum(isnan(group1));
if session==3
    colstodelete=contains(varnames,["DirectionalityIndex","Displacement","Cluster Grade",...
        "DisplacementCorr","Tightness","Incompleteness","StationInTime",...
        "TempMatch","BDistanceClust","BDistanceSpike","nlaps",...
        "rateoverlap","fieldoverlap","lap_perm_stability","stabilityoverlaps",...
        "meanstability","spatialcorrelation"]);
    varnames(colstodelete)=[];
    group1(:,colstodelete)=[];
    group2(:,colstodelete)=[];
elseif session==4
    colstodelete=contains(varnames,["DirectionalityIndex","Cluster Grade",...
        "DisplacementCorr","Tightness","Incompleteness","StationInTime",...
        "TempMatch","BDistanceClust","BDistanceSpike","nlaps",...
        "rateoverlap","fieldoverlap","lap_perm_stability","stabilityoverlaps",...
        "meanstability","spatialcorrelation"]);
    varnames(colstodelete)=[];
    group1(:,colstodelete)=[];
    group2(:,colstodelete)=[];
    
    %     group1(:,contains(varnames,'Displacement'))=abs(wrapTo180(group1(:,contains(varnames,'Displacement'))));
    %     group2(:,contains(varnames,'Displacement'))=abs(wrapTo180(group2(:,contains(varnames,'Displacement'))));
    
end



%% SPLIT BY REGION
% load metadata files and extract region info
group1id=get_region_id(group1id,'D:\Projects\PAE_PlaceCell\AnimalMetadata');
group2id=get_region_id(group2id,'D:\Projects\PAE_PlaceCell\AnimalMetadata');

group1ca1 = group1(strcmp(group1id(:,4),'ca1'),:);
group1ca1id = group1id(strcmp(group1id(:,4),'ca1'),:);
group1ca3 = group1(strcmp(group1id(:,4),'ca3'),:);
group1ca3id = group1id(strcmp(group1id(:,4),'ca3'),:);
group1cortex = group1(strcmp(group1id(:,4),'cortex'),:);
group1cortexid = group1id(strcmp(group1id(:,4),'cortex'),:);

group2ca1 = group2(strcmp(group2id(:,4),'ca1'),:);
group2ca1id = group2id(strcmp(group2id(:,4),'ca1'),:);
group2ca3 = group2(strcmp(group2id(:,4),'ca3'),:);
group2ca3id = group2id(strcmp(group2id(:,4),'ca3'),:);
group2cortex = group2(strcmp(group2id(:,4),'cortex'),:);
group2cortexid = group2id(strcmp(group2id(:,4),'cortex'),:);


[uCA,~,~] = uniqueRowsCA(group1ca1id);
disp([num2str(size(uCA,1)),' control ca1 cells'])
[uCA,~,~] = uniqueRowsCA(group2ca1id);
disp([num2str(size(uCA,1)),' pae ca1 cells'])


[uCA,~,~] = uniqueRowsCA(group1ca3id);
disp([num2str(size(uCA,1)),' control ca3 cells'])
[uCA,~,~] = uniqueRowsCA(group2ca3id);
disp([num2str(size(uCA,1)),' pae ca3 cells'])


[uCA,~,~] = uniqueRowsCA(group1cortexid);
disp([num2str(size(uCA,1)),' control cortex cells'])
[uCA,~,~] = uniqueRowsCA(group2cortexid);
disp([num2str(size(uCA,1)),' pae cortex cells'])


%% place cell filter
%  define_field(group1ca1id,session)

[group1ca1,group1ca1id]=placefieldfilter(group1ca1,group1ca1id,varnames);
[group2ca1,group2ca1id]=placefieldfilter(group2ca1,group2ca1id,varnames);
[group1ca3,group1ca3id]=placefieldfilter(group1ca3,group1ca3id,varnames);
[group2ca3,group2ca3id]=placefieldfilter(group2ca3,group2ca3id,varnames);

[uCA,~,~] = uniqueRowsCA(group1ca1id);
disp([num2str(size(uCA,1)),' control ca1 place cells'])
[uCA,~,~] = uniqueRowsCA(group2ca1id);
disp([num2str(size(uCA,1)),' pae ca1 place cells'])

[uCA,~,~] = uniqueRowsCA(group1ca3id);
disp([num2str(size(uCA,1)),' control ca3 place cells'])
[uCA,~,~] = uniqueRowsCA(group2ca3id);
disp([num2str(size(uCA,1)),' pae ca3 place cells'])

AllStatsca1=CDFplots(group1ca1,group2ca1,{'Sacc','PAE'},varnames,1)
AllStatsca3=CDFplots(group1ca3,group2ca3,{'Sacc','PAE'},varnames,1)


for i=1:length(varnames)
    fig=figure('Name',['ca3 ',varnames{i}],'NumberTitle','off');
    AllStatsca3=CDFplots(group1ca3(:,i),group2ca3(:,i),{'Sacc','PAE'},varnames{i},2);
    toPPT(fig,'exportMode','matlab');
    toPPT('setTitle',AllStatsca3);
    close all
end

visualizecells(uniqueRowsCA(group1ca1id),'control_ca1')
visualizecells(uniqueRowsCA(group2ca1id),'pae_ca1')
visualizecells(uniqueRowsCA(group1ca3id),'control_ca3')
visualizecells(uniqueRowsCA(group2ca3id),'pae_ca3')


%% phase precess
phaseprecess=CDFplots(group1ca1(group1ca1(:,contains(varnames,'Phpval'))<.05,contains(varnames,'PhcircLinCorr')),...
    group2ca1(group2ca1(:,contains(varnames,'Phpval'))<.05,contains(varnames,'PhcircLinCorr')),...
    {'Sacc','PAE'},varnames{contains(varnames,'PhcircLinCorr')},2)

fig=figure;fig.Color=[1 1 1];
[h1,h2]=CoolHistogram(group1ca1(group1ca1(:,contains(varnames,'Phpval'))<.05,contains(varnames,'PhcircLinCorr')),...
    group2ca1(group2ca1(:,contains(varnames,'Phpval'))<.05,contains(varnames,'PhcircLinCorr')),...
    50,varnames{contains(varnames,'PhcircLinCorr')})
xlim([-.4 .4])

print(gcf,'-dpng', '-r400',...
    ['C:\Users\ryanh\Dropbox\school work\UNM\Lab\Projects\PAE Project\Presentations\T32_JournalClub\2019\figures\cylinderPHprec.png'])

group1ca1id(group1ca1(:,contains(varnames,'Phpval'))<.05 & group1ca1(:,contains(varnames,'PhcircLinCorr'))<-.2)

group2ca1id(group2ca1(:,contains(varnames,'Phpval'))<.05 & group2ca1(:,contains(varnames,'PhcircLinCorr'))<-.2)

visualizecells(group1ca1id(group1ca1(:,contains(varnames,'Phpval'))<.05 & group1ca1(:,contains(varnames,'PhcircLinCorr'))<-.01,:),'control')

visualizecells(group2ca1id(group2ca1(:,contains(varnames,'Phpval'))<.05 & group2ca1(:,contains(varnames,'PhcircLinCorr'))<-.1,:),'pae')
%%
% visualizecells(group2ca1id)
% 
% plotexamples(group1ca1id,'control')
% plotexamples(group2ca1id,'pae')
function define_field(groupid,session)
for i=1:length(groupid)
    data=load(groupid{i,1},'spikesID','ratemap','maze_size_cm');
    
    cells=find(contains(data.spikesID.TetrodeNum,groupid(i,2)) & ismember(data.spikesID.CellNum,str2double(groupid(i,3))))';
    upscalefac=15;
    [BW,maskedImage,x,y,fieldarea,X] = segmentImage('map',data.ratemap{cells,session},'figs',true,'upscalefac',upscalefac);
    
%     figure;
%     imagesc(data.ratemap{cells,session});hold on
    for f=1:length(x)
        % downscale
        x_temp=rescale([x{f},upscalefac+1,size(X,2)-upscalefac],1,size(data.ratemap{cells,session},2));
        x{f}=x_temp(1:end-2)';
        y_temp=rescale([y{f},upscalefac+1,size(X,2)-upscalefac],1,size(data.ratemap{cells,session},2));
        y{f}=y_temp(1:end-2)';
        
%         plot(x{f},y{f})

        [majorradius,~,COM] = min_encl_ellipsoid( x{f},y{f});
        
        fields{f}.bounds=[x{f},y{f}];
        fields{f}.peakFR=max(max(data.ratemap{cells,session}(round(y{f}),round(x{f}))));
        fields{f}.COM=COM;
        fields{f}.width=majorradius*2*(data.maze_size_cm(session)/length(data.ratemap{cells,session}));
        fields{f}.fieldarea=fieldarea(f);
    end
    data.openfield=fields;
end
end

function [group,groupid]=placefieldfilter(group,groupid,varnames)
% 1) Minimum peak firing rate of 1 Hz,
% 2) Minimum field width of 8 cm,
% 3) Maximum field width of 80 cm,
% 4) at least 10 trials with consistent behavior.
% 5) at least 100 spikes

groupid=groupid(group(:,contains(varnames,'PeakRate'))>=1 &...
    group(:,contains(varnames,'FieldWidth'))>=9 &...
    group(:,contains(varnames,'FieldWidth'))<=40 &...
    group(:,contains(varnames,'InformationContent'))>=.3 &...
    group(:,contains(varnames,'nSpikes'))>=100,:);

group=group(group(:,contains(varnames,'PeakRate'))>=1 &...
    group(:,contains(varnames,'FieldWidth'))>=9 &...
    group(:,contains(varnames,'FieldWidth'))<=40 &...
    group(:,contains(varnames,'InformationContent'))>=.3 &...
    group(:,contains(varnames,'nSpikes'))>=100,:);
end

function plotexamples(groupid,group)

fig=figure; fig.Color=[1 1 1];
for i=1:size(groupid,1)
    data=load(groupid{i,1});

%         p=postprocessFigures(data,{groupid{i,2},str2double(groupid(i,3))});
    
    cells=find(contains(data.spikesID.TetrodeNum,groupid{i,2}) & ismember(data.spikesID.CellNum,str2double(groupid{i,3})))';
    subplot(3,1,1)
    
    [data_video_spk,~]=createframes_w_spikebinary(data,3-1,cells);
    plot(data_video_spk(:,2),data_video_spk(:,3),'LineWidth',1,'color','k');
    hold on; axis off
    scatter(data_video_spk(data_video_spk(:,6)==1,2),data_video_spk(data_video_spk(:,6)==1,3),10,'r','filled');
    box off; axis image
    title(['nSpikes ',num2str(sum(data_video_spk(:,6)==1))]);
    set(gca,'FontSize',20)
    
    subplot(3,1,2)
    ax = gca;
    SmoothRateMap=data.ratemap{cells,3};
    imAlpha=ones(size(SmoothRateMap));
    imAlpha(isnan(SmoothRateMap))=0;
    imagesc(SmoothRateMap,'AlphaData',imAlpha);
    axis xy; axis off; hold on; box off; axis image;
    colormap(ax,jet(255))
    %                 title([num2str(round(max(max(SmoothRateMap)),2)),' Hz']);
    title(sprintf('IC %4.2f, %1.0fHz',...
        data.measures(cells,contains(data.varnames,["InformationContent","PeakRate"]),3)))
        set(gca,'FontSize',20)

    subplot(3,1,3)
    y=data.thetaautocorr{cells,3};
    plot(y,'LineWidth',2, 'color','k');
    axis tight
    hold on;box off; axis off;axis square
    title(sprintf('Theta %4.1f',...
        data.measures(cells,contains(data.varnames,'thetaindex'),3)))
        set(gca,'FontSize',20)

        set(gcf,'Position',[2 42 273 924])

        print(gcf,'-dpng', '-r400',...
        ['C:\Users\ryanh\Dropbox\school work\UNM\Lab\Projects\PAE Project\Presentations\T32_JournalClub\2019\figures\',group,...
        filesep,groupid{i,1},groupid{i,2},groupid{i,3},'.png'])
    close all
end

end


function visualizecells(groupid,group)
cd('D:\Projects\PAE_PlaceCell\ProcessedData')

sessions=unique(groupid(:,1));

for i=1:length(sessions)
    close all
    idx=find(contains(groupid(:,1),sessions{i}));
    
    data=load(sessions{i});
    
    postprocessFigures.main(data,{groupid(idx,2),str2double(groupid(idx,3))});
    
    
    FigList=findobj(allchild(0), 'flat', 'Type', 'figure');
    for iFig=1:length(FigList)
        FigHandle=FigList(iFig);
        FigName=get(FigHandle, 'Name');
        
        set(FigHandle, 'Position', get(0, 'Screensize'));
        print(gcf,'-dpng', '-r40',...
            ['C:\Users\ryanh\Dropbox\school work\UNM\Lab\Projects\PAE Project\Presentations\SfN2019\CellExamples\cylinder\',group,...
            filesep,groupid{idx(iFig),1},groupid{idx(iFig),2},groupid{idx(iFig),3},'.png'])
    end
    
end
end

