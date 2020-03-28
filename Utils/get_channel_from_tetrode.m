
function channel_num = get_channel_from_tetrode(data,tetrode_path)

% pull out tetrode number to match with channels
[~, filename] = fileparts(tetrode_path);
trodeID = str2double(extractAfter(filename,'TT'));

% find channels
available_channels = data.lfp.channel_list.tetrode_num;

% remove bad channels 
available_channels(~logical(data.lfp.channel_list.good_channels)) = NaN;

% locate first available channel that matches tetrode
channel_num = find(available_channels == trodeID,1,'first');

end