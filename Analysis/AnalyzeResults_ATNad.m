% AnalyzeResults_ATNad LB May 2019
addpath('D:\Users\BClarkLab\GoogleDrive\MatlabDir\BClarkToolbox\Analysis\Visualize')
data=compileResults('F:\ClarkP30_Recordings\ProcessedData');

control={'LB01','LB03','LB05','ATN04','ATN05','ATN08','ATN16','ATN10','ATN17','ATN18'};%,{'LB01','LB03','LB05'}{'ATN04','ATN05','ATN08','ATN16','ATN10'}
transgenic={'LB04','LB06','LB07','ATN07','ATN09','ATN10','ATN14','ATN15'}; %,{'LB02','LB04','LB06','LB07'}{'ATN07','ATN09','ATN14','ATN15'}

%% COMPILE GROUPS
data.control.measures=[]; %control measures
data.control.id=[];
for i=1:length(control)
    data.control.measures=cat(1,data.control.measures,data.(control{i}).measures);
    data.control.id=cat(1,data.control.id,data.(control{i}).id);
end

data.transgenic.measures=[]; %transgenic measures
data.transgenic.id=[];
for i=1:length(transgenic)
    data.transgenic.measures=cat(1,data.transgenic.measures,data.(transgenic{i}).measures);
    data.transgenic.id=cat(1,data.transgenic.id,data.(transgenic{i}).id);
end

% COMPILE DATA
group1=data.control.measures; group2=data.transgenic.measures; group1id=data.control.id; group2id=data.transgenic.id;

%% DELETE MEASURES FOR LINEAR TRACK
varnames=data.varnames; group1(isinf(group1))=NaN; group2(isinf(group2))=NaN;

colstodelete=contains(varnames,["Displacement", "DisplacementCorr","DirectionalityIndex","nlaps","rateoverlap"...
    ,"fieldoverlap","lap_perm_stability","stabilityoverlaps","meanstability","spatialcorrelation"]);

varnames(colstodelete)=[]; group1(:,colstodelete,:)=[]; group2(:,colstodelete,:)=[];

%% HD Filter 
[group1HD,group1HDid]=headdirectionfilter(group1(:,:,:),group1id(:,:,1),varnames);
[group2HD,group2HDid]=headdirectionfilter(group2(:,:,:),group2id(:,:,1),varnames);

%ITERATE THROUGH PROCESSED DATA TO COMPUTE DISTRIBUTIVE RATIO AND
%WITHIN-STABILITY
% DR_Group1=nan(size(group1HDid,1),1,6);
% stability_Group1=nan(size(group1HDid,1),1,6);
% within_Coeff2=nan(size(group2HDid,1),1,6);
% for i=1:size(group1HDid,1)
%     temp=load(['F:\ClarkP30_Recordings\ProcessedData\',group1HDid{i,1}],...
%         'events','frames','spikesID','Spikes','samplerate','ratemap','maze_size_cm');
%     
%     cell=find(contains(temp.spikesID.TetrodeNum,group1HDid{i,2}) & ismember(temp.spikesID.CellNum,str2double(group1HDid{i,3})))';
%     ses=size(temp.events,2);  
%         
%     for ii=1:ses
%     
%         if ses>4
%             continue
%         end
%         
%     [data_video_spk,~]=createframes_w_spikebinary(temp,ii,cell);
%     spks_VEL=data_video_spk(data_video_spk(:,6)==1,4);
%      
%     [~,~,~,~,~,hdTuning]=tuningcurve(data_video_spk(data_video_spk(:,6)==0,4),spks_VEL,temp.samplerate);
%    
%     DR_Group1(i,1,ii) = HD_cell_analysis.distributiveRatio(temp.ratemap{cell,ii},temp.frames,hdTuning,temp.maze_size_cm(ii));
%     
%     stability_Group1(i,1,ii)=HD_cell_analysis.hd_stability(data_video_spk);
%     
%     [within_Coeff1(i,1,ii),~,~] = within_HDstability(data_video_spk,temp.samplerate);
% 
%     end
%     clear temp data_video_spk
% end

% DR_Group2=nan(size(group2HDid,1),1,6);
% stability_Group2=nan(size(group2HDid,1),1,6);
% within_Coeff2=nan(size(group2HDid,1),1,6);
% for i=1:size(group2HDid,1)
%     temp=load(['F:\ClarkP30_Recordings\ProcessedData\',group2HDid{i,1}],...
%         'events','frames','spikesID','Spikes','samplerate','ratemap','maze_size_cm');
%     
%     cell=find(contains(temp.spikesID.TetrodeNum,group2HDid{i,2}) & ismember(temp.spikesID.CellNum,str2double(group2HDid{i,3})))';
%     ses=size(temp.events,2);  
%         
%     for ii=1:ses
%  
%         if ses>4
%             continue
%         end
%     [data_video_spk,~]=createframes_w_spikebinary(temp,ii,cell);
%     spks_VEL=data_video_spk(data_video_spk(:,6)==1,4);
%      
%     [~,~,~,~,~,hdTuning]=tuningcurve(data_video_spk(data_video_spk(:,6)==0,4),spks_VEL,temp.samplerate);
%    
%     DR_Group2(i,1,ii) = HD_cell_analysis.distributiveRatio(temp.ratemap{cell,ii},temp.frames,hdTuning,temp.maze_size_cm(ii));
%     
%     stability_Group2(i,1,ii)=HD_cell_analysis.hd_stability(data_video_spk);
%     
%     [within_Coeff2(i,1,ii),~,~] = within_HDstability(data_video_spk,temp.samplerate);
% 
%     end
%     
%     clear temp
% end
% 
% AllStats=CDFplots(stability_Group1(:,1,1),stability_Group2(:,1,1),{'Control','TgF344-AD'},{'Stability _New'}',1)
% AllStats=CDFplots(within_Coeff1(:,1,1),within_Coeff2(:,1,1),{'Control','TgF344-AD'},{'Stability _Classic'}',1)
% AllStats=CDFplots(DR_Group1(:,1,1),DR_Group2(:,1,1),{'Control','TgF344-AD'},{'Distributive Ratio'}',1)


% AllStats=CDFplots(group1HD(:,:,1),group2HD(:,:,1),{'Control','TgF344-AD'},varnames,1);

%% ___________________________LOCAL FUNCTION BELOW_________________________

function [group,groupid]=headdirectionfilter(group,groupid,varnames)
% 1) Minimum peak firing rate of 1 Hz, 
% 2) Maximum field width of 8 cm,  
% 3) at least 100 spikes
idx=sum(group(:,contains(varnames,'PeakRate'),1)>=1,3)>0 &...
    sum(group(:,contains(varnames,'mean_vector_length'),1)>=.2,3)>0 &...
    sum(group(:,contains(varnames,'nSpikes'),1)>=100,3)>0; % & ...
    %sum(group(:,contains(varnames,'Direct_infoContent'),1)>=.3,3)>0;

groupid=groupid(idx,:);
group=group(idx,:,:);
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

end

end

       