function reorganize_dat_from_chanMap(basepath,channel_map)
% Reorganizes and saves dat file in order of channel map. 

 data.session_path = basepath;
 [~,data.basename] = fileparts(data.session_path);

% Get json path
jsonfile = locate_binary(data.session_path);

% Load mapped data
D = load_open_ephys_binary(jsonfile{1},'continuous',1,'mmap');

% Load xml file 
chanMap = xlsread(channel_map);

% Channel map excel must in in specific format with intan numbers organized
% in column 11. 
channel_num = chanMap(1:64,11)+1; % add 1 for matlab

disp('saving separate .dat files for reorganized headstage data and accelerometer data');
% headstage data
fidout = fopen(fullfile(data.session_path,[data.basename ,'.dat']), 'w');
fwrite(fidout,D.Data.Data.mapped(channel_num,:),'int16');
fclose(fidout);

% Save accelerometer data in separate dat file
data_idx = contains({D.Header.channels.description},'Auxiliar');

% accelerometer data
fidout = fopen(fullfile(data.session_path,[data.basename ,'_accel.dat']), 'w');
fwrite(fidout,D.Data.Data.mapped(data_idx',:),'int16');
fclose(fidout);

end

% locate json file from open ephys output
function lfpfile = locate_binary(session_path)
file = table2cell(struct2table(dir([session_path,filesep,'**\*.oebin'])));
lfpfile = strcat(file(:,2),filesep,file(:,1));
end