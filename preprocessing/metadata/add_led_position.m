% update recording logs with led positions
led_pos = readtable([basepath,filesep,'led_positions.csv']);

files = struct2table(dir([basepath,filesep,'*_metadata.mat']));
for i = 1:length(files.name)
    rat = strsplit(files.name{i},'_');
    if ~ismember(rat{1,1},led_pos.rat_id)
        continue
    end
    idx = find(contains(led_pos.rat_id,rat{1,1}));
    % pull led position 
    led_position.front = led_pos.front(idx);
    led_position.back = led_pos.back(idx);
    led_position.left = led_pos.left(idx);
    led_position.right = led_pos.right(idx);
    %load metadata file
    load([basepath,filesep,rat{1,1},'_metadata.mat']);
    AnimalMetadata.ExtracellEphys.led_position = led_position;
    
    disp('Saving Infomation to .mat...')
    save([basepath,filesep,rat{1,1},'_metadata.mat'],'AnimalMetadata')
    disp('done')
end



