% AnalyzeResults_ATNad LB October 2019
data=compileResults('F:\ClarkP30_Recordings\ProcessedData'); 

% Create groups given rat IDs
control={'LB01','LB03','LB05','ATN04','ATN05','ATN08','ATN14','ATN16','ATN10','ATN17','ATN16','ATN18'};
transgenic={'LB04','LB06','LB07','ATN07','ATN09','ATN15'}; 

% Rats to be included in the analyses based on evidence of electrode
% passing through dorsal hippocampus 
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
met_mat = [];
for x = 1:4
[group_place,groupid_place,place_idx]=place_cell_filter(group,groupid,varnames,x); 
met_mat = [met_mat place_idx];
end

%Creating indicies to use later
group_place = group(met_mat(:,1) == 1,:,:); groupid_place = groupid(met_mat(:,1) == 1,:); place_critera = met_mat(met_mat(:,1)== 1,:);
genotype = double(contains(groupid_place(:,1),transgenic)); % 1 if Tg+, 0 if control; 
region = double(contains(groupid_place(:,1),'ATN')); % 1 if ATN, 0 if cortex; 

% Quick stats run through
stat_plot(group_place(genotype==0,:,1),group_place(genotype==1,:,1),{'WT','Tg+'},varnames)

% Save "place cell" figures to folder to view & check 
visualize_cells(groupid(place_idx,:),'d:\Users\BClarkLab\Google Drive (lberkowitz@unm.edu)\Projects\ClarkP30_Ephys\Analysis\Place_data\Figures\Met_Place')

% compile measures for place cells only 
all_place_data = [region genotype ones(size(genotype,1),1) group_place(:,:,1);...
    region genotype ones(size(genotype,1),1)+1 group_place(:,:,2);...
    region genotype ones(size(genotype,1),1)+2 group_place(:,:,3);...
    region genotype ones(size(genotype,1),1)+3 group_place(:,:,4)]; 

% Convert cell to table to prep for exporting to csv
varnames = [{'region','genotype','sess_num'} varnames]; varnames = regexprep(varnames, '\W', '');
place_data_all = cell2table(num2cell(all_place_data),'VariableNames',varnames);

% Save data to csv to use in R 
writetable(place_data_all,'d:\Users\BClarkLab\Google Drive (lberkowitz@unm.edu)\Projects\ClarkP30_Ephys\Analysis\Place_data\place_data_all.csv')


ses = 1; %set relative to the first session 
for i = 1:length(groupid_place)
    
    temp=load(['F:\ClarkP30_Recordings\ProcessedData\',groupid_place{i,1}],...
        'events','frames','spikesID','Spikes','samplerate','ratemap','maze_size_cm');
    
    cell=find(contains(temp.spikesID.TetrodeNum,groupid_place{i,2}) & ismember(temp.spikesID.CellNum,str2double(groupid_place{i,3})))';
    
    %Stability measure of how similar each ratemap is between the first two
    %standard sessions 
    if size(temp.events,2) == 4 && place_critera(i,ses) == 1 && place_critera(i,ses+1) == 1
        RateMap1 = temp.ratemap{cell,ses};
        RateMap2 = temp.ratemap{cell,ses+1};
        RateMap1(isnan(RateMap1))=0; RateMap1(isinf(RateMap1))=0;
        RateMap1 = padarray(RateMap1,[3 3],0,'both');
        RateMap2(isnan(RateMap2))=0; RateMap2(isinf(RateMap2))=0;
        RateMap2 = padarray(RateMap2,[3 3],0,'both');
        similarity(i,1) = corr2(RateMap1,RateMap2);
    else
        similarity(i,1) = nan; 
    end
    
    % Lets take a look at rotation 
    if size(temp.events,2) == 4 && place_critera(i,ses) == 1 && place_critera(i,ses+1) == 1
        RateMap1 = temp.ratemap{cell,ses};
        RateMap2 = temp.ratemap{cell,ses+1};
        [Rotation_mat(i,1),Rotation_cor_mat(i,1)] = Displacement2(RateMap1,RateMap2);
    else
        Rotation_mat(i,1) = nan; Rotation_cor_mat(i,1) = nan;
    end
    
    
    clear RateMap1 RateMap2 temp cell
%     fig=figure('Name',id{i,1},'NumberTitle','off');
%     postprocessFigures.ratemaps_2d(fig,temp.ratemap{cell,ses},group_place(i,:,ses),varnames)
%    
%     export_fig('d:\Users\BClarkLab\Google Drive (lberkowitz@unm.edu)\Projects\ClarkP30_Ephys\Analysis\Place_data\Figures','-eps'
end

%Histogram of stability scores for place cells found with 4 condition
%recordings
figure; 
stat_plot(similarity(genotype==0,:),similarity(genotype==1,:),{'WT','Tg+'},{'Spatial Correlation'},'plots',2)

fig = figure;
fig.Color = [1 1 1];
subplot(1,2,1)
h1=polarhistogram(deg2rad(Rotation_mat(genotype == 0,1)),'BinWidth',.1047)
h1.FaceColor = [.5 .5 .5]
title('Degree of Peak Ratemap Correlation, F344')
subplot(1,2,2)
h=polarhistogram(deg2rad(Rotation_mat(genotype == 1,1)),'BinWidth',.1047)
h.FaceColor = 'r'
title('Degree of Peak Ratemap Correlation, TgF344-AD')

%Run a v-test to determine if shift angles are significantly clustered
%around predicted angle (0 or 90 deg for stability and rotation,
%respectively).
[WT_pval WT_v] = circ_vtest(deg2rad(Rotation_mat(genotype == 0 & ~isnan(Rotation_mat(:,1)),1)), deg2rad(90))

%Run a watson U to determine if the two distributions are significantly
%overlapping. 
[p,U2_obs,U2_H0]=watsons_U2_perm_test(deg2rad(Rotation_mat(genotype == 0 & ~isnan(Rotation_mat(:,1)),1)),deg2rad(Rotation_mat(genotype == 1 & ~isnan(Rotation_mat(:,1)),1)),100)

%% ___________________________LOCAL FUNCTION BELOW_________________________


function [group,groupid,idx]=place_cell_filter(group,groupid,varnames,ses)
% 1) Minimum peak firing rate of 1 Hz,
% 2) Maximum field width of 8 cm,
% 3) at least 100 spikes
idx=group(:,contains(varnames,'InformationContent'),ses) >= .3 & ...
    group(:,contains(varnames,'FieldWidth'),ses) >= 9 & ...
    group(:,contains(varnames,'FieldWidth'),ses) <= 40 & ...
    group(:,contains(varnames,'PeakRate'),ses) >= 1 & ...
    group(:,contains(varnames,'nSpikes'),ses) >= 100;

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


