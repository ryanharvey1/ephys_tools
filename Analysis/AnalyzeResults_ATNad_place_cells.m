% AnalyzeResults_ATNad LB October 2019
addpath('D:\Users\BClarkLab\ephys_tools\external_packages\mvmdist-master')
data=compileResults('F:\ClarkP30_Recordings\ProcessedData');

control={'LB01','LB03','LB05','ATN04','ATN05','ATN08','ATN14','ATN16','ATN10','ATN17','ATN16','ATN18'};
transgenic={'LB04','LB06','LB07','ATN07','ATN09','ATN15'}; 

rats={'ATN05','ATN07','ATN08','ATN09','ATN10','ATN14','ATN15','ATN17','ATN16','ATN17','ATN18'};

%% COMPILE DATA FROM INDIVIDUAL RATS
data.measures=[]; %control measures
data.id=[];
for i=1:length(rats)
    data.measures=cat(1,data.measures,data.(rats{1,i}).measures);
    data.id=cat(1,data.id,data.(rats{i}).id);
end

%% COMPILE DATA
group=data.measures; groupid=data.id;

%% DELETE MEASURES FOR LINEAR TRACK
varnames=data.varnames; group(isinf(group))=NaN;
colstodelete=contains(varnames,["Displacement", "DisplacementCorr","DirectionalityIndex","nlaps","rateoverlap"...
    ,"fieldoverlap","lap_perm_stability","stabilityoverlaps","meanstability","spatialcorrelation"]);
varnames(colstodelete)=[]; group(:,colstodelete,:)=[];

clear colstodelete data

%% Forgo analyzing cells with low firing rate (< 1hz) and limited spikes (< 100spikes) in the first session
[group,groupid]=quality_filter(group,groupid,varnames);

%% Compile data from above with Measures 
group(:,:,5:end) = []; %remove data from exploratory sessions

%% place cell filter 
[group_place,groupid_place,place_idx]=place_cell_filter(group,groupid,varnames); 
genotype = double(contains(groupid_place(:,1),transgenic)); % 1 if Tg+, 0 if control; 
region = double(contains(groupid_place(:,1),'ATN')); % 1 if ATN, 0 if cortex; 

%Lets look at what the place cell filter finds
data = load('ATN07_S20180712182441','openfield','ratemap')
figure;
nm = ceil(sqrt(size(data.ratemap,1)));
for i = 1:size(data.ratemap,1)
    subplot(nm,nm,i)
    imAlpha=ones(size(data.ratemap{i,1}));
    imAlpha(isnan(data.ratemap{i,1}))=0;
    imagesc(data.ratemap{i,1},'AlphaData',imAlpha);
    axis xy; axis off; hold on; box off; axis image;
    colormap(viridis(255))
    data.openfield{1, 1}.fields{1, i}
    hold on
    for f = 1:length(data.openfield{1, 1}.fields{1, i}.bounds)
        plot(data.openfield{1, 1}.fields{1, i}.bounds{f}(:,1),...
            data.openfield{1, 1}.fields{1, i}.bounds{f}(:,2),'LineWidth',3)
    end
    pause(.00001)
end

stat_plot(group_place(genotype==0,:,1),group_place(genotype==1,:,1),{'WT','Tg+'},varnames)

% visualize_cells(groupid(place_idx,:),'d:\Users\BClarkLab\Google Drive (lberkowitz@unm.edu)\Projects\ClarkP30_Ephys\Analysis\Place_data\Figures\Met_Place')

% compile measures
all_data = [region genotype ones(size(genotype,1),1)  group_place(:,:,1);...
    region genotype ones(size(genotype,1),1)+1  group_place(:,:,2);...
    region genotype ones(size(genotype,1),1)+2 group_place(:,:,3);...
    region genotype ones(size(genotype,1),1)+3 group_place(:,:,4)]; 

varnames = [{'region','genotype','sess_num'} varnames];
place_data_all = cell2table(num2cell(all_data),'VariableNames',varnames);

% writetable(place_data_all,'d:\Users\BClarkLab\Google Drive (lberkowitz@unm.edu)\Projects\ClarkP30_Ephys\Analysis\Place_data\place_data_all.csv')

%% Some example figs 
% plot heatplots for top 20 IC per group

%Gather data for each group
IC_Tg = place_data_all(genotype == 1,contains(varnames,'corrected_info_content'));
IC_WT = place_data_all(genotype == 0,contains(varnames,'corrected_info_content'));

%Find top 20 location
[~,Tg_idx]=sort(IC_Tg.InformationContent,1);
[~,WT_idx]=sort(IC_WT.InformationContent,1);


Tg_place_id = groupid_place(Tg_idx <= 20,:);
ses = 1; %only look at session 1 
for i = 1:length(Tg_place_id)
    
    temp=load(['F:\ClarkP30_Recordings\ProcessedData\',Tg_place_id{i,1}],...
        'events','frames','spikesID','Spikes','samplerate','ratemap','maze_size_cm');
    
    cell=find(contains(temp.spikesID.TetrodeNum,Tg_place_id{i,2}) & ismember(temp.spikesID.CellNum,str2double(Tg_place_id{i,3})))';
        
    [data_video_spk,~]=createframes_w_spikebinary(temp,ses,cell);
    [~,within,~]=HD_cell_analysis.four_quarter_stability(data_video_spk,temp.samplerate,'std');
        
    % Plot ratemap 
    
    fig = figure;
    postprocessFigures.ratemaps_2d(fig,temp.ratemap{i,ses},group_place(i,:,ses),varnames)
    xlabel(Tg_place_id{i,1})
%     
%     print(fig,'-dpng', '-r300',['d:\Users\BClarkLab\Desktop\Laura Temp\TG_HD',...
%         HD_cell_id{i,1},'.png'])
%     close
%     
    
end

% Phase precession (spike on phase) plots for Tg+ and WT rats 


%% ___________________________LOCAL FUNCTION BELOW_________________________

function [group,groupid,idx]=headdirectionfilter(group,groupid,varnames)
% 1) Minimum peak firing rate of 1 Hz,
% 2) Maximum field width of 8 cm,
% 3) at least 100 spikes
idx=group(:,contains(varnames,'Direct_infoContent'),1) >= .2 & ...
    group(:,contains(varnames,'nSpikes'),1) >= 100 & ...
    group(:,contains(varnames,'mean_vector_length'),1) >= .2 & ...
    group(:,contains(varnames,'PeakRate'),1) >= 1;

groupid=groupid(idx,:);
group=group(idx,:,:);
end

function [group,groupid,idx]=bidirectional_filter(group,groupid,varnames)
% 1) Minimum peak firing rate of 1 Hz,
% 2) Maximum field width of 8 cm,
% 3) at least 100 spikes
idx=group(:,contains(varnames,'num_components'),1) == 2 & ...
    group(:,contains(varnames,'Sparsity'),1) >= .3 & ...
    group(:,contains(varnames,'PeakRate'),1) >= 1 & ...
    group(:,contains(varnames,'nSpikes'),1) >= 100;

groupid=groupid(idx,:);
group=group(idx,:,:);
end

function [group,groupid,idx]=place_cell_filter(group,groupid,varnames)
% 1) Minimum peak firing rate of 1 Hz,
% 2) Maximum field width of 8 cm,
% 3) at least 100 spikes
idx=group(:,contains(varnames,'InformationContent'),1) >= .3 & ...
    group(:,contains(varnames,'FieldWidth'),1) >= 9 & ...
    group(:,contains(varnames,'FieldWidth'),1) <= 40 & ...
    group(:,contains(varnames,'PeakRate'),1) >= 1 & ...
    group(:,contains(varnames,'nSpikes'),1) >= 100;

groupid=groupid(idx,:);
group=group(idx,:,:);
end

function [group,groupid]=quality_filter(group,groupid,varnames)
% 1) Minimum peak firing rate of 1 Hz,
% 2) at least 100 spikes

idx=group(:,contains(varnames,'PeakRate'),1) >= 1 &...
    group(:,contains(varnames,'nSpikes'),1) >= 100;

groupid=groupid(idx,:);
group=group(idx,:,:);
end

function [group,groupid,idx]=interneuron_filter(group,groupid,varnames)
% 1) spike width less than .2  
% 2) overall firing rate more than 10z

idx=group(:,contains(varnames,'OverallFiringRate'),1) >= 10 & ... 
    group(:,contains(varnames,'spikewidth'),1) <= .2 ;

groupid=groupid(idx,:);
group=group(idx,:,:);
end

%% ADD TO INDICATED SPOT AT LATER TIME AFTER EFFICIENCY IS IMPROVED (CURRENT CODE TAKES TOO MUCH MEMORY TO RUN - LB 10/11/2019

%         num_spikes = group(i,contains(varnames,'nSpikes'),ii);
                
%         % only fit vMMM for predicted components to be less than 4 & has at least 100 spikes.
%         if num_components == 2 && num_spikes > 100 
            
%             for model = 1:2
%                 % get tuning curve to create samples
%                 
%                 % obtain the raw tuning curve
%                 tuning = temp.hdTuning{cell,ii};
%                 
%                 % round and multiply to obtain an estimate of degrees per bin
%                 dir_prop = round(tuning*10);
%                 
%                 % create a distribution of angles that reflects the tuning
%                 % curve proportion.
%                 samples_deg=[];
%                 for d = 1:size(bin_centers,2)
%                     samples_deg = [samples_deg; repmat(bin_centers(1,d),dir_prop(1,d),1)];
%                 end
%                 
%                 % fit angular distribution with vmm
%                 samples=wrapToPi(deg2rad(samples_deg));
%                 fittedVmm = fitmvmdist(samples, model,'MaxIter',300, 'Replicates', 10);
%                 
%                 %obtain the fitted angles
%                 angles = linspace(-pi, pi, size(samples,1))';
%                 likelihoodsFitted = fittedVmm.pdf(angles);
%                 
%                 % obtain the mean squared error to estimate the error between
%                 % fitted and raw angular distribution
%                 err(model)=immse(likelihoodsFitted,samples);
%                 
%                 % store the loglikelihood for each model
%                 fit_measures(model)=fittedVmm.logLikelihood;
%                 
%                 
%             end
%             
%             % compare loglikelihoods - must have improvement of 20% over
%             % next simplest model to win. If that is not met, winner is
%             % most simplest model (unimodal)
%             winner = find([diff(fit_measures),0]./abs(fit_measures)>=.2);
%             if isempty(winner)
%                 winner = 1;
%             end
%             vmm_winner = winner(1);
%             
%             % Lets see if the autocorr method (from get_nComponents) is in
%             % agreement with the winner.
%             if vmm_winner == num_components
%                 winner = num_components;
%                 loglike = fit_measures(winner);
%                 err = err(winner);
%                 fit = 1; 
%             else
%                 winner = NaN;
%                 loglike = NaN;
%                 err = NaN;
%                 fit = 1;
%             end
%             
%             clear samples_deg
%         else 
%                 winner = NaN;
%                 loglike = NaN;
%                 err = NaN;
%                 fit = 0;
%         end
%         
%         
%         hd_data(i,6,ii) = winner;
%         hd_data(i,9,ii) = loglike;
%         hd_data(i,8,ii) = err;
%         hd_data(i,7,ii) = fit;
        % Lets calculate stability
%         [data_video_spk,~]=createframes_w_spikebinary(temp,ii,cell);
%         
%         stability=HD_cell_analysis.stability(data_video_spk,temp.samplerate);
%         
%         hd_data(i,6,ii) = stability;
        

% TgHD_1 = first(contains(HD_cell_id(:,1),transgenic),:);
% [~,pdIdx]=max(TgHD_1,[],2);
%     [~,I]=sort(pdIdx);
%     TgHD_1=TgHD_1(I,:);
%     TgHD_2=second(contains(HD_cell_id(:,1),transgenic),:);
%     TgHD_3=third(contains(HD_cell_id(:,1),transgenic),:);
%     TgHD_4=fourth(contains(HD_cell_id(:,1),transgenic),:);
%     
%     figure; 
%     subaxis(1,4,1)
%     imagesc(TgHD_1)
%     colormap(viridis(255))
%     subaxis(1,4,2)
%     imagesc(TgHD_2)
%     colormap(viridis(255))
%     subaxis(1,4,3)
%     imagesc(TgHD_3)
%     colormap(viridis(255))
%     subaxis(1,4,4)
%     imagesc(TgHD_4)
%     colormap(viridis(255))
% 
% 
% WtHD_1 = first(~contains(HD_cell_id(:,1),transgenic),:);
% [~,pdIdx]=max(WtHD_1,[],2);
%     [~,I]=sort(pdIdx);
%     WtHD_1=WtHD_1(I,:);
%     WtHD_2=second(~contains(HD_cell_id(:,1),transgenic),:);
%     WtHD_3=third(~contains(HD_cell_id(:,1),transgenic),:);
%     WtHD_4=fourth(~contains(HD_cell_id(:,1),transgenic),:);
%     
%     figure; 
%     subaxis(1,4,1)
%     imagesc(WtHD_1)
%     colormap(viridis(255))
%     subaxis(1,4,2)
%     imagesc(WtHD_2)
%     colormap(viridis(255))
%     subaxis(1,4,3)
%     imagesc(WtHD_3)
%     colormap(viridis(255))
%     subaxis(1,4,4)
%     imagesc(WtHD_4)
%     colormap(viridis(255))