% Probe pipeline 
addpath(genpath('D:\Users\BClarkLab\ephys_tools\external_packages\analysis-tools'))


% CD to your data 

%.0  Reorganize dat file if your data was not mapped when streamed to disk
reorganize_dat_from_chanMap(pwd,'CED_E1_4X16_front_front.xlsx')

% 1. Get filtered raw dat signal
filter_raw_dat_from_dat

% 2. Get lfp, initial xml, and channel map. 
get_lfp

% 3. In Neuroscope, check channel map, skip bad channels, and save. 

% 4. Update channel map from .xml
create_channelmap_assy156(pwd)

% 5. Copy data over to SSD 
[~,basename] = fileparts(pwd);
command = ['robocopy ',pwd,' ','C:\Kilo_temp',filesep,basename];
system(command);

% 5. Spike sort using kilosort 2.0
run_ks2(['C:\Kilo_temp',filesep,basename])

% 6. Copy back over to data file 
disp('Saving data back to data folder from ssd')
command = ['robocopy ','C:\Kilo_temp',filesep,basename,' ',pwd,' /e'];
system(command);

% 7. Spike sort in Phy

% 4. After spike sort cleanup 
after_spikesort_cleanup.handle_phy

% 5. Posprocess 

