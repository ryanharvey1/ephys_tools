% AnalyzeResults_ATNad LB October 2019
addpath('D:\Users\BClarkLab\ephys_tools\external_packages\mvmdist-master')
data=compileResults('F:\ClarkP30_Recordings\ProcessedData');

control={'LB03','LB05','ATN04','ATN05','ATN08','ATN14','ATN16','ATN10','ATN17','ATN16','ATN18'};
transgenic={'LB04','LB06','LB07','LB10','ATN07','ATN09','ATN15'}; 

rats={'LB03','LB04','LB05','LB06','LB07','ATN04','ATN05','ATN07','ATN08','ATN09','ATN10','ATN14','ATN15','ATN17','ATN16','ATN18'};

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

% Initalize Matrix for new measures
hd_data_vars = {'cell_id','session_num','num_components','mean_prominence','score','four_quarter','signal2noise'};
hd_data = zeros(size(groupid,1),size(groupid,2),size(groupid,3));

% LOOP THROUGH CELLS TO FIND ESTIMATED NUMBER OF PEAKS 
for i=1:size(groupid,1)
    temp=load(['F:\ClarkP30_Recordings\ProcessedData\',groupid{i,1}],...
        'events','frames','spikesID','Spikes','samplerate','ratemap','hdTuning');
    
    cell=find(contains(temp.spikesID.TetrodeNum,groupid{i,2}) & ismember(temp.spikesID.CellNum,str2double(groupid{i,3})))';
    
    ses=size(temp.events,2);
   
    
    for ii = 1:ses
        
        if ii > 4 %don't look at the exploratory data sessions
            continue
        end
        
        hd_data(i,1,ii) = cell;
        hd_data(i,2,ii) = ii;
        
        % Let's check for multimodal cells using an autocor methods adapted
        % from Kevin Allen's group. 
        [num_components , mean_prom , score] = get_nComponents(temp,ii, cell,0);
        
        hd_data(i,3,ii) = num_components;
        hd_data(i,4,ii) = mean_prom;
        hd_data(i,5,ii) = max(score);
   
        [data_video_spk,~]=createframes_w_spikebinary(temp,ii,cell);
        [within_Coeff,~,~]=HD_cell_analysis.four_quarter_stability(data_video_spk,temp.samplerate,'std');
        
        hd_data(i,6,ii)  = within_Coeff;
        
        pfr = group(cell,contains(varnames, 'PeakRate'),ii);
        
        afr = group(cell,contains(varnames, 'OverallFiringRate'),ii);
        
        signal2noise = pfr./afr;
        
        hd_data(i,7,ii) = signal2noise; 
       
        clear data_video_spk 
        

    end

    
    
end

%% Compile data from above with Measures 
group(:,:,5:6) = []; %remove data from exploratory sessions

% Cat HD data with postprocess measures 
dir_measures_all = cat(2,hd_data,group);
vars = [hd_data_vars varnames];


fig = figure; 
fig.Color = [1 1 1];
r = dir_measures_all(:,contains(vars, 'mean_vector_length'),1);
dic = dir_measures_all(:,contains(vars, 'Direct_infoContent'),1);
peakfr = dir_measures_all(:,contains(vars, 'PeakRate'),1);

scatter(r,dic,'filled','r')
hold on; 
ax=gca;
plot([.2 .2],[ax.YLim(1) ax.YLim(2)],'--k','LineWidth',2)
plot([ax.XLim(1) ax.XLim(2)],[.2 .2],'--k','LineWidth',2)

ylabel('Direct Info Content (Bits/Spike)')
xlabel('mean vector length')
set(gca,'FontSize',14,'FontName','Helvetica')


%% HD Filter
[groupHD,groupHDid,HD_idx]=headdirectionfilter(dir_measures_all,groupid,vars);
% visualize_cells(groupHDid,'d:\Users\BClarkLab\Desktop\Laura Temp\HD_met')

%% Bidirectional filter 
[group_BD,groupid_BD,BD_idx]=bidirectional_filter(dir_measures_all,groupid,vars);

%% place cell filter 
[group_place,groupid_place,place_idx]=place_cell_filter(dir_measures_all,groupid,vars);

%% Interneuron filter 

 [group,groupid_IN,IN_idx]=interneuron_filter(dir_measures_all,groupid,vars);

 IN_idx_final = IN_idx & ~HD_idx & ~BD_idx & ~place_idx; %can't meet critera for spatial cells. 
 HD_idx_final = HD_idx & ~BD_idx & ~place_idx & ~IN_idx; 
 BD_idx_final = BD_idx & ~HD_idx & ~place_idx & ~IN_idx; 
 place_idx_final = place_idx & ~HD_idx & ~BD_idx & ~IN_idx; 
 
 visualize_cells(groupid(BD_idx_final,:),'d:\Users\BClarkLab\Desktop\Laura Temp\BD_met')
 
% Create region ID 
region = double(contains(groupid(:,1),'ATN')); % 1 if ATN, 0 if Cortex 
genotype = double(contains(groupid(:,1),transgenic)); % 1 if Tg+, 0 if control; 
vars = [{'region','genotype','HD_cell','BD_cell','place_Cell','IN_cell'} vars];

% compile measures
all_data = [region genotype HD_idx_final BD_idx_final place_idx_final IN_idx_final dir_measures_all(:,:,1);...
    region genotype HD_idx_final BD_idx_final place_idx_final IN_idx_final dir_measures_all(:,:,2);...
    region genotype HD_idx_final BD_idx_final place_idx_final IN_idx_final dir_measures_all(:,:,3);...
    region genotype HD_idx_final BD_idx_final place_idx_final IN_idx_final dir_measures_all(:,:,4)]; 

dir_id = [groupid;groupid;groupid;groupid];

HD_data_all = cell2table([dir_id num2cell(all_data)],'VariableNames',[{'ID','tetrode','channel'} vars]);

writetable(HD_data_all,'D:\Users\BClarkLab\Google Drive (lberkowitz@unm.edu)\UNM PhD\Conferences\SFN 2019\data\hd_data_all.csv')

%% Some example figs 

HD_cell_data = dir_measures_all(HD_idx & ~place_idx & ~BD_idx & ~genotype,:);
HD_cell_id = groupid(HD_idx & ~place_idx & ~BD_idx & ~genotype,:);
for i = 1:length(HD_cell_id)
    
    temp=load(['F:\ClarkP30_Recordings\ProcessedData\',HD_cell_id{i,1}],...
        'events','frames','spikesID','Spikes','samplerate','ratemap','maze_size_cm','hdTuning');
    
    cell=find(contains(temp.spikesID.TetrodeNum,HD_cell_id{i,2}) & ismember(temp.spikesID.CellNum,str2double(HD_cell_id{i,3})))';
    
    ses=size(temp.events,2);
    
    [data_video_spk,~]=createframes_w_spikebinary(temp,1,cell);
    [~,within,~]=HD_cell_analysis.four_quarter_stability(data_video_spk,temp.samplerate,'std');
        
    % Plot Spikes on Path by HD 
    fig = figure;
    
    fig.Color = [1 1 1];
    
    angBins=0:6:360;
    bin_centers=movmedian(angBins,2);
    bin_centers(1)=[];
    Polarplot = polar(deg2rad(bin_centers),temp.tuning{cell,1},'b');
    set(Polarplot,'linewidth',1,'color','k');
    axis off
    set(0,'Showhiddenhandles','on')
    extrastuff = setdiff(get(gca,'children'),Polarplot);
    delete(extrastuff)
    hold on
    horizontal=line([-max(tuning) max(tuning)],[0 0]); % for running max and min
    vertical=line([0 0],[-max(tuning) max(tuning)]);
    set(horizontal,'linewidth',2,'color',[.4 .4 .4]);
    set(vertical,'linewidth',2,'color',[.4 .4 .4]);
    axis image
    
    uistack(horizontal,'bottom')
    uistack(vertical,'bottom')
    
    title(sprintf('r: %4.2f DIC: %4.2f' ,[r,Ispk]))
    
    h=fill(get(Polarplot, 'XData'), get(Polarplot, 'YData'),...
        'k');
    set(h,'FaceAlpha',.5)
    
    postprocessFigures.plot_HD_tuning(temp,1,cell)
    xlabel(HD_cell_id{i,1})
%     
%     print(fig,'-dpng', '-r300',['d:\Users\BClarkLab\Desktop\Laura Temp\TG_HD',...
%         HD_cell_id{i,1},'.png'])
%     close
%     
    
end


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