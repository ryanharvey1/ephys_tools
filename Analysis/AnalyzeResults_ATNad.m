% AnalyzeResults_ATNad

addpath('D:\Users\BClarkLab\GoogleDrive\MatlabDir\BClarkToolbox\Analysis\Visualize')
data=compileResults('F:\ClarkP30_Recordings\ProcessedData');

control={'LB03','LB05'};%,'LB01','LB03','LB05'};'ATN04','ATN05',,'LB05''LB03''ATN08'
transgenic={'LB04','LB06'}; %,'LB02','LB04','LB06'};'ATN01','ATN03','LB04','LB06','ATN07','ATN09','ATN10','ATN14'

%% COMPILE GROUPS
data.control.measures=[];
data.control.id=[];
for i=1:length(control)
    data.control.measures=cat(1,data.control.measures,data.(control{i}).measures);
    data.control.id=cat(1,data.control.id,data.(control{i}).id);
end

data.transgenic.measures=[];
data.transgenic.id=[];
for i=1:length(transgenic)
    data.transgenic.measures=cat(1,data.transgenic.measures,data.(transgenic{i}).measures);
    data.transgenic.id=cat(1,data.transgenic.id,data.(transgenic{i}).id);
end

% COMPILE DATA
group1=data.control.measures;
group2=data.transgenic.measures;
group1id=data.control.id;
group2id=data.transgenic.id;


%% DELETE MEASURES FOR LINEAR TRACK
varnames=data.varnames;
group1(isinf(group1))=NaN;
group2(isinf(group2))=NaN;

colstodelete=contains(varnames,["Displacement", "DisplacementCorr","DirectionalityIndex","nlaps","rateoverlap"...
    ,"fieldoverlap","lap_perm_stability","stabilityoverlaps","meanstability","spatialcorrelation"]);

varnames(colstodelete)=[];
group1(:,colstodelete,:)=[];
group2(:,colstodelete,:)=[];


%% HD Filter 
[group1HD,group1HDid]=headdirectionfilter(group1(:,:,:),group1id(:,:,1),varnames);
[group2HD,group2HDid]=headdirectionfilter(group2(:,:,:),group2id(:,:,1),varnames);

AllStats=CDFplots(group1HD(:,:,1),group2HD(:,:,1),{'Control','TgF344-AD'},varnames,1);

tuning1=get_tuning(group1HDid);
tuning2=get_tuning(group2HDid);





%Get first quarter index for Controls
[~,pdIdx]=max(tuning1.S1,[],2);
[~,ATNidx]=sort(pdIdx);
ATN_sorted=tuning1.S1(ATNidx,:);
ATN_sorted_2=tuning1.S2(ATNidx,:);
ATN_sorted_3=tuning1.S3(ATNidx,:);
ATN_sorted_4=tuning1.S4(ATNidx,:);



%Get first quarter index for Transgenics
[~,pdIdx]=max(tuning2.S1,[],2);
[~,ATNidx]=sort(pdIdx);
ATN_sorted2=tuning2.S1(ATNidx,:);
ATN_sorted2_2=tuning2.S2(ATNidx,:);
ATN_sorted2_3=tuning2.S3(ATNidx,:);
ATN_sorted2_4=tuning2.S4(ATNidx,:);

vars=fieldnames(tuning1);
for i=1:size(vars,1)
    temp=tuning1.(vars{i});
    for ii=1:size(temp,2)
    cells(ii,:,i)=rescale(temp,ii,ii+1)
    end
end




%Create PopVec Figure
fig=figure;
fig.Color=[1 1 1];
subaxis(2,4,1)
imagesc(ATN_sorted);axis off
title('Session 1')
subaxis(2,4,2)
imagesc(ATN_sorted_2);axis off
title('Session 2')
subaxis(2,4,3)
imagesc(ATN_sorted_3);axis off
title('Session 3 - Rotated')
subaxis(2,4,4)
imagesc(ATN_sorted_4);axis off
title('Session 4')

subaxis(2,4,5)
imagesc(ATN_sorted2);axis off
subaxis(2,4,6)
imagesc(ATN_sorted2_2);axis off
subaxis(2,4,7)
imagesc(ATN_sorted2_3);axis off
subaxis(2,4,8)
imagesc(ATN_sorted2_4); axis off

y=[group1HD(:,9,1) group1HD(:,9,2) group1HD(:,9,3) group1HD(:,9,4)]; x=(1:size(y,2));
SEM=nanstd(y)/sqrt(size(y,1));
hold on; 
shadedErrorBar(x,nanmean(y),SEM,'-k',1); 

y=[group2HD(:,9,1) group2HD(:,9,2) group2HD(:,9,3) group2HD(:,9,4)]; x=(1:size(y,2));
SEM=nanstd(y)/sqrt(size(y,1));
shadedErrorBar(x,nanmean(y),SEM,'-r',1); 
% ylim([.3 1]);

data=[group1HD(:,9,1); group1HD(:,9,2); group1HD(:,9,3); group1HD(:,9,4); nan(1,1)];
distIdx=[ones(size(group1HD,1),1); ones(size(group1HD,1),1); ones(size(group1HD,1),1);...
    ones(size(group1HD,1),1); 2];
catIdx=[ones(size(group1HD,1),1); ones(size(group1HD,1),1).*2; ones(size(group1HD,1),1).*3;...
    ones(size(group1HD,1),1).*4;1];

plotSpread(data,'categoryIdx',distIdx,'distributionIdx',catIdx,...
    'categoryMarkers',{'o','o'},'categoryColors',{[.25 .25 .25],'k'})

data=[group2HD(:,9,1); group2HD(:,9,2); group2HD(:,9,3); group2HD(:,9,4);nan(1,1)];
distIdx=[ones(size(group2HD,1),1);ones(size(group2HD,1),1) ...
    ;ones(size(group2HD,1),1);ones(size(group2HD,1),1);2];
catIdx=[ones(size(group2HD,1),1);ones(size(group2HD,1),1).*2 ...
    ;ones(size(group2HD,1),1).*3;ones(size(group2HD,1),1).*4;1];

plotSpread(data,'categoryIdx',distIdx,'distributionIdx',catIdx,...
    'categoryMarkers',{'o','o'},'categoryColors',{'r','k'})


%All Data - R length Beeswarm with clusters that have min 100 spikes
data=[group1(group1(:,39)>100,9,1); group1(group1(:,39)>100,9,2);...
    group1(group1(:,39)>100,9,3); group1(group1(:,39)>100,9,4); nan(1,1)];
distIdx=[ones(size(group1(group1(:,39)>100,9,1),1),1); ...
    ones(size(group1(group1(:,39)>100,9,1),1),1); ...
    ones(size(group1(group1(:,39)>100,9,1),1),1);...
    ones(size(group1(group1(:,39)>100,9,1),1),1); 2];
catIdx=[ones(size(group1(group1(:,39)>100,9,1),1),1); ...
    ones(size(group1(group1(:,39)>100,9,1),1),1).*2; ...
    ones(size(group1(group1(:,39)>100,9,1),1),1).*3;...
    ones(size(group1(group1(:,39)>100,9,1),1),1).*4;1];

plotSpread(data,'categoryIdx',distIdx,'distributionIdx',catIdx,...
    'categoryMarkers',{'o','o'},'categoryColors',{[.25 .25 .25],'k'})

data=[group2(group2(:,39)>100,9,1); group2(group2(:,39)>100,9,2);...
    group2(group2(:,39)>100,9,3); group2(group2(:,39)>100,9,4);nan(1,1)];
distIdx=[ones(size(group2(group2(:,39)>100,9,1),1),1);...
    ones(size(group2(group2(:,39)>100,9,1),1),1) ...
    ;ones(size(group2(group2(:,39)>100,9,1),1),1);...
    ones(size(group2(group2(:,39)>100,9,1),1),1);2];
catIdx=[ones(size(group2(group2(:,39)>100,9,1),1),1);...
    ones(size(group2(group2(:,39)>100,9,1),1),1).*2 ...
    ;ones(size(group2(group2(:,39)>100,9,1),1),1).*3;...
    ones(size(group2(group2(:,39)>100,9,1),1),1).*4;1];

plotSpread(data,'categoryIdx',distIdx,'distributionIdx',catIdx,...
    'categoryMarkers',{'o','o'},'categoryColors',{'r','k'})

fig=figure;
figure.Color=[1 1 1];
subaxis(1,4,1)
[h1,h2]=CoolHistogram(group1(group1(:,39)>100,9,1),group2(group2(:,39)>100,9,1),10,'mean Vector Length')

subaxis(1,4,2)
[h1,h2]=CoolHistogram(group1(group1(:,39)>100,9,2),group2(group2(:,39)>100,9,2),10,'mean Vector Length')

subaxis(1,4,3)
[h1,h2]=CoolHistogram(group1(group1(:,39)>100,9,3),group2(group2(:,39)>100,9,3),10,'mean Vector Length')

subaxis(1,4,4)
[h1,h2]=CoolHistogram(group1(group1(:,39)>100,9,4),group2(group2(:,39)>100,9,4),10,'mean Vector Length')

%% ___________________________LOCAL FUNCTION BELOW_________________________
function [group,groupid]=placefieldfilter(group,groupid,varnames)
% 1) Minimum peak firing rate of 2 Hz, 
% 2) Minimum field width of 8 cm, 
% 3) Maximum field width of 180 cm, 
% 4) Mean vector length less than .5 to limit direction x place. 
% 5) at least 100 spikes
%IDX HD cells
% idx=sum(group(:,contains(varnames,'PeakRate'),1)>=1,3)>0 &...
%     sum(group(:,contains(varnames,'mean_vector_length'),1)>=.2,3)>0 &...
%     sum(group(:,contains(varnames,'nSpikes'),1)>=100,3)>0;

%IDX place cells
idx=sum(group(:,contains(varnames,'PeakRate'),1)>=2,3)>0 &...
    sum(group(:,contains(varnames,'FieldWidth'),1)>=8,3)>0 &...
    sum(group(:,contains(varnames,'FieldWidth'),1)<=40,3)>0 & ...
    sum(group(:,contains(varnames,'InformationContent'),1)>=.4,3)>0 & ...
    sum(group(:,contains(varnames,'nSpikes'),1)>=100,3)>0;

groupid=groupid(idx,:);
group=group(idx,:,:);

end

%%
function [group,groupid]=headdirectionfilter(group,groupid,varnames)
% 1) Minimum peak firing rate of 1 Hz, 
% 2) Maximum field width of 8 cm,  
% 3) at least 100 spikes
idx=sum(group(:,contains(varnames,'PeakRate'),1)>=1,3)>0 &...
    sum(group(:,contains(varnames,'mean_vector_length'),1)>=.3,3)>0 &...
    sum(group(:,contains(varnames,'nSpikes'),1)>=100,3)>0; % & ...
    %sum(group(:,contains(varnames,'Direct_infoContent'),1)>=.3,3)>0;

groupid=groupid(idx,:);
group=group(idx,:,:);
end


function [group,groupid]=qualityfilter(group,groupid,varnames)
% 1) Minimum peak firing rate of 2 Hz, 
% 2) Maximum field width of 8 cm,  
% 3) at least 100 spikes
idx=sum(group(:,contains(varnames,'ShortISI'),1)<.05,3)>0 &...
    sum(group(:,contains(varnames,'nSpikes'),1)>=100,3)>0;

groupid=groupid(idx,:);
group=group(idx,:,:);
end

function [group,groupid]=interneuronfilter(group,groupid,varnames)
% Only keep putative interneurons 
%
% 1.peak to valley duration above .2ms
% 2.thetaindex about 2
idx=sum(group(:,contains(varnames,'spikewidth'),1)<.2,3)>0 &...
    sum(group(:,contains(varnames,'thetaindex'),1)>=2,3)>0 &...
    sum(group(:,contains(varnames,'nSpikes'),1)>200,3)>0;

groupid=groupid(idx,:);
group=group(idx,:,:);
end

function visualizecells(groupid,group)
cd('F:\ClarkP30_Recordings\ProcessedData')
for i=1:length(groupid)
    data=load(groupid{i,1});
    
    p=postprocessFigures(data,{groupid{i,2},str2double(groupid(i,3))});
    set(gcf, 'Position', get(0, 'Screensize'));
    
    print(gcf,'-dpng', '-r300',...
        ['d:\Users\BClarkLab\Desktop\Laura Temp\PlaceCells\',group,...
        filesep,groupid{i,1},groupid{i,2},groupid{i,3},'.png'])
    close all
    
end
end


function autoCorrMat=popVecAutocorr(groupid)
cd('F:\ClarkP30_Recordings\ProcessedData')
autoCorrMat=[];
for i=1:length(groupid)
    data=load(groupid{i,1},'thetaautocorr','spikesID');
    %     close all
    cells=find(contains(data.spikesID.TetrodeNum,groupid{i,2}) & ismember(data.spikesID.CellNum,str2double(groupid(i,3))))';
    
    for ii=cells
        autoCorrMatTEMP=data.thetaautocorr{cells,1}';
        autoCorrMat=[autoCorrMat; autoCorrMatTEMP];  
    end
end

end

function tuning=get_tuning(groupid)
cd('F:\ClarkP30_Recordings\ProcessedData')

sessions={'S1','S2','S3','S4'};
for i=1:length(groupid)
    data=load(groupid{i,1},'events','spikesID','Spikes','frames','samplerate');
    cells=find(contains(data.spikesID.TetrodeNum,groupid{i,2}) & ismember(data.spikesID.CellNum,str2double(groupid(i,3))))';
    
    if size(data.events,2)<4
        continue
    end
    % GET TUNNING CURVE
    for ii=cells
        for iii=1:4
            [data_video_spk,~]=createframes_w_spikebinary(data,iii,ii);
            [~,~,~,~,~,temptuning]=tuningcurve(data_video_spk(data_video_spk(:,6)==0,4),...
                data_video_spk(data_video_spk(:,6)==1,4),data.samplerate);
            if exist('tuning','var') && isfield(tuning,(sessions{iii})) 
                r=size(tuning.(sessions{iii}),1)+1;
            else
                r=1;
            end
            tuning.(sessions{iii})(r,:)=rescale(temptuning,0,1);
        end 
    end
   
 
%     polarplot(deg2rad(0:6:360),tuning,'k')
%     ax=gca;
%     set(ax,'RGrid','off','ThetaGrid','off','ThetaTick',[0 90],...
%         'ThetaTickLabels',[0 90],'RTick',max(tuning),'RTickLabel',max(tuning),'LineWidth',3);
%     axis tight
end

end

% FOR DIAG POP VECTOR
function [mat]=arrangenorm(mat)
fields=fieldnames(mat);
[~,I]=max(mat.(fields{1}),[],2);
[~,I2]=sort(I);
for i=1:length(fields)
    mat.(fields{i})=mat.(fields{i})(I2,:);
end
end
       