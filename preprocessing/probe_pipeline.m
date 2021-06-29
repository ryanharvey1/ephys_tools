 Probe pipeline 
%
% Procedures for processing electrophysiology data using Cambridge Neurotech
% probes assy 156 in the BClark lab. 
%
% Dependencies
%   - ephys_tools>external_packages>buzcode 
%
% LBerkowitz 2020

% Set paths for data 
data_folder = pwd;
[~,basename] = fileparts(data_folder);
probe_map = 'CED_E1_4X16_front_front.xlsx';
ssd_path = 'C:\Kilo_temp';

addpath(genpath('D:\Users\BClarkLab\ephys_tools\external_packages\analysis-tools'))

%% ##########################################################################

%                       Setup data for clustering  

% ##########################################################################

%% Preprocessing (organize dat based on channel map) 
% Channels were not organized during acquisition. This opens and saves the
% dat file so channels are organized based on channel map. 

%.0  Reorganize dat file if your data was not mapped when streamed to disk
reorganize_dat_from_chanMap(data_folder,probe_map)

% 1. Get filtered raw dat signal
filter_raw_dat_from_dat

% 2. Apply CAR to remove hypersynchronous events 
applyCARtoDat(['filtered','.dat'], 64, data_folder);

% 3. Make xml file
write_xml(data_folder,'channel_map_64.xlsx',30000,1250) 

%% 4. In Neuroscope, check channel map, skip bad channels, and save. 

% follow steps for chosen spike sorting method (Kilosort2.5 or Klusta)

%% ##########################################################################

%                       Sorting using Kilosort2.5

% ##########################################################################

%% For Kilosort: Create lfp file and make channel map
% 5. Get lfp, initial xml, and channel map. 
get_lfp

% 6. Update channel map from basename.xml
create_channelmap_assy156(data_folder)

%% Spike sorting
% Move chanMap, xml, and dat to SSD folder. 
%creating a folder on the ssd for chanmap,dat, and xml
ks2_folder = fullfile(ssd_path, basename);
mkdir(ks2_folder);

% 7. Spike sort using kilosort 2.5 (data on ssd)
run_ks2([ssd_path,filesep,basename])

% 8. Copy back over to data file 
disp('Saving data back to data folder from ssd')
command = ['robocopy ',ssd_path,filesep,basename,' ',data_folder,' /e'];
system(command);

%% 9. Clean up kilo results in Phy

%% 10. After spike sort cleanup 
after_spikesort_cleanup.handle_phy

% 11. Posprocess 
postprocess('DLC',1)

%% ##########################################################################

%                       Sorting using Klusta 

% ##########################################################################
%% For Klusta
d   = dir([data_folder,filesep,'*.xml']);
parameters = LoadXml(fullfile(data_folder,d(1).name));

% 5. create klusta_folders as well as prb and prm files 
makeProbeMap(data_folder)

% 6. Runs klusta on each shank by calling system
run_klusta(data_folder)

%% 7. Spike sort in Phy (Spike sorting kwik files in phy will update the kwik file)

%% 9. After spike sort cleanup 
after_spikesort_cleanup.handle_kwik(data_folder)

% 10. Posprocess 
postprocess('DLC',1)
