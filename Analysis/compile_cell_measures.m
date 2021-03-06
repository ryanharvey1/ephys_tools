% compile_cell_measures LB August 2020
% Compiles all measures across all sessions from indicated animals (rats).
% Cells with less that 100 spikes and < 1Hz peak rate are excluded. Saves
% compiled data to all_cells.csv and groupid for cells as cell_list.csv stored in location indicated by the save_path
% variable. 
function compile_cell_measures(rats,transgenic)
% Input:
%  rats: cell array of animal IDs to include e.g. {'LB03','LB10'}
if ~exist('save_path','var')
    save_path = 'F:\ClarkP30_Recordings\Analysis\';
end

%% Pull Data from ProcessedData
data=compileResults('F:\ClarkP30_Recordings\ProcessedData');

%% Compile Data for individual rats
data.measures=[]; %control measures
data.id=[];
for i=1:length(rats)
    data.measures=cat(1,data.measures,data.(rats{1,i}).measures);
    data.id=cat(1,data.id,data.(rats{i}).id);
end

%% Unpack structure into matrix and cell array, respectively.
group=data.measures; groupid=data.id;

%% Delete Measures for linear Track
varnames=data.varnames; group(isinf(group))=NaN;
colstodelete=contains(varnames,["Displacement", "DisplacementCorr","DirectionalityIndex","nlaps","rateoverlap"...
    ,"fieldoverlap","lap_perm_stability","stabilityoverlaps","meanstability","spatialcorrelation"]);
varnames(colstodelete)=[]; group(:,colstodelete,:)=[];

% Address gaps in variable names that prevent concatenation when making
% 'vars' below
varnames = regexprep(varnames, ' ', '');

clear colstodelete data

%% Forgo analyzing cells with low firing rate (< 1hz) and limited spikes (< 100spikes) in the first session
[group,groupid]=quality_filter(group,groupid,varnames);

cell_list = table('Size',[size(groupid)],...
    'VariableTypes',{'cell','cell','double'},'VariableNames',{'SessionID','tetrode','cell'});
cell_list.SessionID = groupid(:,1); cell_list.tetrode = groupid(:,2); cell_list.cell = groupid(:,3);

writetable(cell_list,[save_path,'cell_list.csv'])

%% Create region ID 
region = double(contains(groupid(:,1),'ATN')); % 1 if ATN, 0 if Cortex 
genotype = double(contains(groupid(:,1),transgenic)); % 1 if Tg+, 0 if control; 
vars = [{'region','genotype','Condition_num'} varnames];

%% Compile 
all_data = [region genotype ones(size(group,1),1) group(:,:,1);...
    region genotype ones(size(group,1),1)+1 group(:,:,2);...
    region genotype ones(size(group,1),1)+2 group(:,:,3);...
    region genotype ones(size(group,1),1)+3 group(:,:,4)]; 

dir_id = [groupid;groupid;groupid;groupid];

%% Save Data
data_all = cell2table([dir_id num2cell(all_data)],'VariableNames',[{'SessionID','tetrode','channel'} vars]);
writetable(data_all,[save_path,'all_cells.csv'])

clearvars -except data_all cell_list save_path

end
%% ___________________________LOCAL FUNCTION BELOW_________________________


function [group,groupid]=quality_filter(group,groupid,varnames)
% 1) Minimum peak firing rate of 1 Hz,
% 2) at least 100 spikes

idx=group(:,contains(varnames,'PeakRate'),1) >= 1 &...
    group(:,contains(varnames,'nSpikes'),1) >= 100;

groupid=groupid(idx,:);
group=group(idx,:,:);
end

%% ADD TO INDICATED SPOT AT LATER TIME AFTER EFFICIENCY IS IMPROVED (CURRENT CODE TAKES TOO MUCH MEMORY TO RUN - LB 10/11/2019
% % Initalize Matrix for new measures
% hd_data_vars = {'cell_id','session_num','num_components','mean_prominence','score','four_quarter','signal2noise'};
% hd_data = zeros(size(groupid,1),size(groupid,2),size(groupid,3));
% 
% % LOOP THROUGH CELLS TO FIND ESTIMATED NUMBER OF PEAKS 
% for i=1:size(groupid,1)
%     temp=load(['F:\ClarkP30_Recordings\ProcessedData\',groupid{i,1}],...
%         'events','frames','spikesID','Spikes','samplerate','ratemap','hdTuning');
%     
%     cell=find(contains(temp.spikesID.TetrodeNum,groupid{i,2}) & ismember(temp.spikesID.CellNum,str2double(groupid{i,3})))';
%     
%     ses=size(temp.events,2);
%    
%     
%     for ii = 1:ses
%         
%         if ii > 4 %don't look at the exploratory data sessions
%             continue
%         end
%         
%         hd_data(i,1,ii) = cell;
%         hd_data(i,2,ii) = ii;
%         
%         % Let's check for multimodal cells using an autocor methods adapted
%         % from Kevin Allen's group. 
%         [num_components , mean_prom , score] = get_nComponents(temp,ii, cell,0);
%         
%         hd_data(i,3,ii) = num_components;
%         hd_data(i,4,ii) = mean_prom;
%         hd_data(i,5,ii) = max(score);
%    
%         [data_video_spk,~]=createframes_w_spikebinary(temp,ii,cell);
%         [within_Coeff,~,~]=HD_cell_analysis.four_quarter_stability(data_video_spk,temp.samplerate,'std');
%         
%         hd_data(i,6,ii)  = within_Coeff;
%         
%         pfr = group(cell,contains(varnames, 'PeakRate'),ii);
%         
%         afr = group(cell,contains(varnames, 'OverallFiringRate'),ii);
%         
%         signal2noise = pfr./afr;
%         
%         hd_data(i,7,ii) = signal2noise; 
%        
%         clear data_video_spk 
%         
% 
%     end
% 
%     
%     
% end

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