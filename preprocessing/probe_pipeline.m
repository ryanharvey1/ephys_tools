% Probe pipeline 
addpath(genpath('D:\Users\BClarkLab\ephys_tools\external_packages\analysis-tools'))

% CD to your data 
[~,basename] = fileparts(pwd);
%.0  Reorganize dat file if your data was not mapped when streamed to disk
reorganize_dat_from_chanMap(pwd,'CED_E1_4X16_front_front.xlsx')

% 1. Apply CAR to remove high frequency noise 
applyCARtoDat([basename ,'.dat'], 64, pwd);

% 2. Get filtered raw dat signal
filter_raw_dat_from_dat

% 3. Make xml file
write_xml(pwd,'channel_map_64.xlsx',30000,1250) 
% 4. In Neuroscope, check channel map, skip bad channels, and save. 

% 5. Get lfp, initial xml, and channel map. 
get_lfp

% 6. In Neuroscope, check channel map, skip bad channels, and save. 

% 7. Update channel map from .xml
create_channelmap_assy156(pwd)

% 8. Copy data over to SSD 
command = ['robocopy ',pwd,' ','C:\Kilo_temp',filesep,basename];
system(command);

% 9. Spike sort using kilosort 2.0
run_ks2(['C:\Kilo_temp',filesep,basename])

% 10. Copy back over to data file 
disp('Saving data back to data folder from ssd')
command = ['robocopy ','C:\Kilo_temp',filesep,basename,' ',pwd,' /e'];
system(command);

% 11. Spike sort in Phy

% 12. After spike sort cleanup 
after_spikesort_cleanup.handle_phy

% 13. Posprocess 
postprocess('DLC',1)
