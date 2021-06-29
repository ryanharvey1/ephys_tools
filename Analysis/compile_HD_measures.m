% Compile_HD_measures_ATNad.m: compiles measures computed in postprocess.m
% and computes additional measures assessing HD cell function.


% Set Groups
control = {'LB03','PoS4','LB11'};
transgenic = {'LB04','LB07','LB10','PoS3'};

% Set Rats to include in analysis
rats = {'LB03','LB04','LB07','LB10','LB11','PoS3','PoS4'};

%% Load a list of classified HD cells for analysis
% Grab list of classified cells
groupid = readtable('d:\Users\BClarkLab\github\Berkowitz_et_al_2021\Cell_Classification\data\cell_list.csv');
baseline_index = groupid.baseline_index;

% Initalize variables
vars = {'cell_idx','Condition_num','PeakRate',...
    'nSpikes','r','dic','stability',...
    'sig2noise','pref_dir','kappa','r_test_p','r_test_z'};

% Initalize measure matrix
group = zeros(size(groupid,1),size(vars,2)); % 8 is equal to number of sessions

% Bin Centers
da=pi/30; 
bin_centers=da/2:da:2*pi-da/2;

%% Loop through cells and compute HD measures
for i = 1 : size(groupid.session,1)
    
    % Load data file and find < *Data*Session*Cell >
    data = load(['F:\ClarkP30_Recordings\ProcessedData\',groupid.session{i,1}],...
        'events','frames','spikesID','Spikes','samplerate','ratemap','hdTuning',...
        'hdTuning_corrected','measures','varnames');
    
    % Pull things we need from data structure;
    samplerate = data.samplerate;
    varnames = data.varnames;
    cell = groupid.cell_idx(i);
    
    % Loop over conditions (changed to just look at baseline lb 3/30/21)
    %     for ii = baseline_index(i)
    
    ii = 1;
    
    measures = data.measures(:,:,baseline_index(i));
    hdTuning = data.hdTuning{cell,baseline_index(i)};
    
    % Get Frames
    [data_video_spk,~]=createframes_w_spikebinary(data,baseline_index(i),cell);
    
    [tuning_overtime,stability] = HD_cell_analysis.stability(data_video_spk,samplerate);
    
    
    
    % Build group with measures
    group(i,1,ii)  =  cell;                                                      % cell idx
    group(i,2,ii)  =  baseline_index(i);                                                        % condition number
    group(i,3,ii)  =  measures(cell,contains(varnames,'PeakRate'));              % Peak Firing Rate
    group(i,4,ii)  =  measures(cell,contains(varnames,'nSpikes'));               % tuning stable
    group(i,5,ii)  =  measures(cell,contains(varnames,'mean_vector_length'));    % mean vector length
    group(i,6,ii)  =  measures(cell,contains(varnames,'Direct_infoContent'));     % directional info content
    group(i,7,ii)  =  stability;                                                 % tuning stable
    
    % Smooth Tuning Curve
    if ~isnan(hdTuning)
        smooth_tune = smooth_hd([hdTuning hdTuning hdTuning]); %smooth over 36 degrees
        smooth_tune = smooth_tune(1, 61:120);
        
        % Signal to noise Peak Rate/Min Rate
        [peak,ang] = max(smooth_tune); % grab
        group(i,8,ii) =  peak/mean(smooth_tune);
        
        % Pull Peak Direction
        group(i,9,ii) = rad2deg(bin_centers(ang));
        
        % compute kappa
        group(i,10,ii) = circ_kappa(bin_centers',hdTuning');
        
        % perform Rayleigh test for circular distribution
        [group(i,11,ii), group(i,12,ii)] = circ_rtest(bin_centers', smooth_tune);
        
    else
        group(i,8,ii) = NaN;
        group(i,9,ii) = NaN;
        group(i,10,ii) = NaN;
        group(i,11,ii) = NaN;
        group(i,12,ii) = NaN;
    end
    
    clear data_video_spk measures
    
    
    %     end
end

%% Get Condition Type from metadata
rats = extractBefore(groupid.session,'_');
sess = extractBefore(extractAfter(groupid.session,'_'),'.mat');

% Loop through and categorize
for i = 1:length(rats)
    
    load([rats{i},'_metadata.mat'])
    sess_temp = sess{i};
    
    if AnimalMetadata.RecordingLogs.(sess_temp).MazeTypes == "Cylinder" || ...
            AnimalMetadata.RecordingLogs.(sess_temp).MazeTypes == "pedestal,Cylinder,pedestal"
        condition(i,1) = {'baseline'};
        
    elseif AnimalMetadata.RecordingLogs.(sess_temp).MazeTypes == "Cylinder,Cylinder"
        condition(i,1) = {'stability'};
        
    elseif AnimalMetadata.RecordingLogs.(sess_temp).MazeTypes == "Cylinder,Cylinder,Cylinder"
        condition(i,1) = {'dim'};
        
    elseif AnimalMetadata.RecordingLogs.(sess_temp).MazeTypes == "Cylinder,Cylinder,Cylinder,Cylinder" || ...
            AnimalMetadata.RecordingLogs.(sess_temp).MazeTypes == "Cylinder,Cylinder,Cylinder,Cylinder,box" || ...
            AnimalMetadata.RecordingLogs.(sess_temp).MazeTypes == "Cylinder,Cylinder,Cylinder,Cylinder,circ track" || ...
            AnimalMetadata.RecordingLogs.(sess_temp).MazeTypes == "Cylinder,Cylinder,Cylinder,Cylinder,Cylinder" || ...
            AnimalMetadata.RecordingLogs.(sess_temp).MazeTypes == "pedestal,Cylinder,Cylinder,Cylinder,Cylinder,pedestal" || ...
            AnimalMetadata.RecordingLogs.(sess_temp).MazeTypes == "Cylinder,Cylinder,Cylinder,Cylinder,Cylinder,Cylinder"
        condition(i,1) = {'stability_rotation'};
        
    elseif AnimalMetadata.RecordingLogs.(sess_temp).MazeTypes == "Cylinder,box,Cylinder" || ...
            AnimalMetadata.RecordingLogs.(sess_temp).MazeTypes == "Cylinder,box" || ...
            AnimalMetadata.RecordingLogs.(sess_temp).MazeTypes == "pedestal,Cylinder,box,Cylinder,pedestal"
        condition(i,1) = {'box'};
        
    end
end

%% Add Region, Genotype
% Create region ID
region = double(contains(groupid.session,'ATN')); % 1 if ATN, 0 if Cortex
genotype = double(contains(groupid.session,transgenic)); % 1 if Tg+, 0 if control;
vars = [{'session','Condition_type','region','genotype'} vars];

% %% Arrange into long format data frame
% all_data = [region genotype group(:,:,1);...
%     region genotype  group(:,:,2);...
%     region genotype group(:,:,3);...
%     region genotype group(:,:,4)];
% id = [groupid(:,1), condition ;groupid(:,1), condition ;groupid(:,1), condition ;groupid(:,1), condition ];
%% Arrange into long format data frame for baseline only
all_data = [region genotype group(:,:,1)];

id = [groupid.session, condition ];

%% Save Data as CSV
HD_data_all = cell2table([id, num2cell(all_data)],'VariableNames',vars);

writetable(HD_data_all,'d:/Users/BClarkLab/github/Berkowitz_et_al_2021/HD_Modulation/hd_metrics_baseline.csv')


function smoothed = smooth_hd(tuning)
gaus = @(tuning,mu,sig)exp(-(((tuning-mu).^2)/(2*sig.^2)));
window = 10;
x = -window/2:window/2;
mu = 0;
sig = 3;
kernel = gaus(x,mu,sig);
smoothed = conv(tuning,kernel/sum(kernel),'same');
end




