% get_area_per_channel

data_path = 'F:\Projects\PAE_PlaceCell\ProcessedData\';

sessions = dir(fullfile(data_path,'*.mat'));
sessions = struct2table(sessions);

group_id = [];
for i = 1:length(sessions.name)
data = load(fullfile(data_path,sessions.name{i}),'session_path');
 
session_info = LoadParameters(data.session_path);

r=length(session_info.channels);
channels = cellstr([repmat('TT',r,1),num2str(session_info.channels'+1),repmat('.mat',r,1)]);
channels = regexprep(channels, '\s', '');

temp_id = [cellstr(repmat(sessions.name{i},r,1)),...
    channels,...
    channels];

group_id = [group_id;temp_id];
end

group_id=get_region_id(group_id,'F:\Projects\PAE_PlaceCell\AnimalMetadata');
% make dg into ca3
for i=1:length(group_id)
    if strcmp(group_id{i,end},'dg')
        group_id{i,end} = 'ca3';
    end
end

df = cell2table(group_id,'VariableNames',{'session' 'channel' 'cell' 'area'});
df.channel = str2double(extractBetween(df.channel,'TT','.mat'));
df.session = extractBefore(df.session,'.mat');
df.cell = [];

writetable(df,'F:\Projects\PAE_PlaceCell\analysis\channel_list.csv')
