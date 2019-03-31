% 
% Created by Ben C June 2014

% identify path to data files
path = '/Users/bjclark/Desktop/Dropbox/ec recordings/conversion/bc217/2010-07-19_16-55-26';
ReadData = FindFiles('bc217_TT4_u1_s04.txt', 'StartingDirectory', path);

for i = 1:length(ReadData);

    data = importdata(ReadData{i}); % load xy data

    % extract coords, spikes, angle, and direction from ReadData output
    datai = data.data(:,2); % red x-coords
    datai(:,2) = data.data(:,3); % red y-coords
    datai(:,3) = data.data(:,4); % green x-coords
    datai(:,4) = data.data(:,5); % green y-coords

    % find red LED non-detects
    NDs = find(datai(:,1) == 255 & datai(:,2) == 255);
    NDsFILT = datai(NDs,:);
    numbNDs = numel(NDsFILT);

    keep NDs NDsFILT ReadData path i numbNDs
    [filepath, filename] = fileparts(ReadData{i});
    save([filepath filesep filename '_NDs.mat']);
    
end