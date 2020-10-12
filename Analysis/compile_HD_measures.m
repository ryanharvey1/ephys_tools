% Compile_HD_measures_ATNad.m: compiles measures computed in postprocess.m
% and computes additional measures assessing HD cell function.
addpath('F:\ClarkP30_Recordings\AnimalMetaData');
addpath('F:\ClarkP30_Recordings\Analysis\')

% Set Groups
control = {'LB03'};
transgenic = {'LB04','LB07','LB10'};

% Set Rats to include in analysis
rats = {'LB03','LB04','LB07','LB10'};

%% Load a list of classified HD cells for analysis
% Grab list of classified cells
list = table2cell(readtable('F:\ClarkP30_Recordings\Analysis\hd_cell_list.csv'));

% set as groupid to use in loop below
groupid = list;

% Initalize variables
vars = {'cell_idx','Condition_num','PeakRate',...
    'nSpikes','r','dic','stability',...
    'sig2noise','pref_dir'};

% Initalize measure matrix
group = zeros(size(groupid,1),size(vars,2),6); % 6 is equal to number of sessions

% Bin Centers
bin_centers = movmedian(0:6:360,2);

%% Loop through cells and compute HD measures
for i = 1 : size(groupid,1)
    
    % Load data file and find < *Data*Session*Cell >
    data = load(['F:\ClarkP30_Recordings\ProcessedData\',groupid{i,1}],...
        'events','frames','spikesID','Spikes','samplerate','ratemap','hdTuning',...
        'hdTuning_corrected','measures','varnames');
    ses = size(data.events,2);
    tet = sscanf(groupid{i,2},'TT%d.mat');
    cell = find_cells(data,tet,groupid{i,3});
    
    % Pull things we need from data structure;
    samplerate = data.samplerate;
    varnames = data.varnames;
    
    
    % Loop over conditions
    for ii = 1 : ses
        
        measures = data.measures(:,:,ii);
        hdTuning = data.hdTuning_corrected{cell,ii}';
        if ii > 4 %don't look at the exploratory data sessions
            continue
        end
        
        % Get Frames
        [data_video_spk,~]=createframes_w_spikebinary(data,ii,cell);

        stability = HD_cell_analysis.stability(data_video_spk,samplerate);
        
        % Build group with measures
        group(i,1,ii)  =  cell;                                                      % cell idx
        group(i,2,ii)  =  ii;                                                        % condition number
        group(i,3,ii)  =  measures(cell,contains(varnames,'PeakRate'));              % Peak Firing Rate
        group(i,4,ii)  =  measures(cell,contains(varnames,'nSpikes'));               % tuning stable
        group(i,5,ii)  =  measures(cell,contains(varnames,'mean_vector_length'));    % mean vector length
        group(i,6,ii)  =  measures(cell,contains(varnames,'Direct_infoContent'));     % directional info content
        group(i,7,ii)  =  stability;                                                 % tuning stable
        
        % Smooth Tuning Curve
        if ~isnan(hdTuning)
            smooth_tune = smoothdata([hdTuning hdTuning hdTuning],'gaussian',6); %smooth over 36 degrees
            smooth_tune = smooth_tune(1, 61:120);
            
            % Signal to noise Peak Rate/Min Rate
            [peak,ang] = max(smooth_tune); % grab
            group(i,8,ii) =  peak/min(smooth_tune);                                    % signal to noise
            
            % Pull Peak Direction
            group(i,9,ii) = bin_centers(ang); % multiply by bin size
        else
            group(i,8,ii) = NaN;                                    % signal to noise
            group(i,9,ii) = NaN; % multiply by bin size
        end
        clear data_video_spk measures
        
        
    end
end

%% Get Condition Type from metadata
rats = extractBefore(groupid(:,1),'_');
sess = extractBefore(extractAfter(groupid(:,1),'_'),'.mat');

% Loop through and categorize
for i = 1:length(rats)
    
    load([rats{i},'_metadata.mat'])
    sess_temp = sess{i};
    
    if AnimalMetadata.RecordingLogs.(sess_temp).MazeTypes == "Cylinder"
        condition(i,1) = {'baseline'};
        
    elseif AnimalMetadata.RecordingLogs.(sess_temp).MazeTypes == "Cylinder,Cylinder"
        condition(i,1) = {'stability'};
        
    elseif AnimalMetadata.RecordingLogs.(sess_temp).MazeTypes == "Cylinder,Cylinder,Cylinder"
        condition(i,1) = {'dim'};
        
    elseif AnimalMetadata.RecordingLogs.(sess_temp).MazeTypes == "Cylinder,Cylinder,Cylinder,Cylinder" || ...
            AnimalMetadata.RecordingLogs.(sess_temp).MazeTypes == "Cylinder,Cylinder,Cylinder,Cylinder,box" || ...
            AnimalMetadata.RecordingLogs.(sess_temp).MazeTypes == "Cylinder,Cylinder,Cylinder,Cylinder,circ track" || ...
            AnimalMetadata.RecordingLogs.(sess_temp).MazeTypes == "Cylinder,Cylinder,Cylinder,Cylinder,Cylinder"
        condition(i,1) = {'stability_rotation'};
        
    elseif AnimalMetadata.RecordingLogs.(sess_temp).MazeTypes == "Cylinder,box,Cylinder"
        condition(i,1) = {'delta_context'};
        
    end
end

%% Add Region, Genotype
% Create region ID
region = double(contains(groupid(:,1),'ATN')); % 1 if ATN, 0 if Cortex
genotype = double(contains(groupid(:,1),transgenic)); % 1 if Tg+, 0 if control;
vars = [{'SessionID','Condition_type','region','genotype'} vars];

%% Arrange into long format data frame
all_data = [region genotype group(:,:,1);...
    region genotype  group(:,:,2);...
    region genotype group(:,:,3);...
    region genotype group(:,:,4)];

id = [groupid(:,1), condition ;groupid(:,1), condition ;groupid(:,1), condition ;groupid(:,1), condition ];
%% Save Data as CSV
HD_data_all = cell2table([id, num2cell(all_data)],'VariableNames',vars);

writetable(HD_data_all,'F:\ClarkP30_Recordings\Analysis\hd_data_all.csv')


