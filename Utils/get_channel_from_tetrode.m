
function channel_num = get_channel_from_tetrode(data,cell)

% pull out tetrode number to match with channels
trodeID = str2double(extractBetween(data.spikesID.TetrodeNum{cell},'TT','.mat'));

% find channels
available_channels = data.lfp.channel_list.tetrode_num;

% remove bad channels 
available_channels(~logical(data.lfp.channel_list.good_channels)) = NaN;

% locate first available channel that matches tetrode
channel_num = find(available_channels == trodeID,1,'first');

end