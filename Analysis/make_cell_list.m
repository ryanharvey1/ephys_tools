% make_cell_list

data_path = 'F:\ClarkP30_Recordings\ProcessedData\';

sessions = dir(fullfile(data_path,'*.mat'));
sessions = struct2table(sessions);
group_id = [];
for i = 1:length(sessions.name)
    data = load(sessions.name{i},'spikesID');

    r=length(data.spikesID.TetrodeNum);

    temp_id = [cellstr(repmat(sessions.name{i},r,1)),...
        data.spikesID.TetrodeNum,...
        cellstr(num2str(data.spikesID.CellNum))];
    group_id = [group_id;temp_id];
end
group_id=get_region_id(group_id,'F:\ClarkP30_Recordings\AnimalMetaData');
% make dg into ca3
for i=1:length(group_id)
    if strcmp(group_id{i,end},'dg')
        group_id{i,end} = 'ca3';
    end
end
group_id = cell2table(group_id,'VariableNames',{'session' 'tt' 'cell' 'area'});

group_id.tt = str2double(extractBetween(group_id.tt,'TT','.mat'));
group_id.session = extractBefore(group_id.session,'.mat');
group_id.cell = str2double(group_id.cell);


writetable(group_id,'F:\ClarkP30_Recordings\analysis/cell_list.csv')
