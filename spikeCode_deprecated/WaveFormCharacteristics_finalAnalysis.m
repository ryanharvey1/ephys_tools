

path = '/Users/bjclark/Documents/MATLAB/test_data/spike_data/waveforms';

f_list = dir([path '/*.mat']);

for i = 1:length(f_list);
    meanDur_all(i) = load(f_list(i).name, 'meanDur');
    maxPeak_all(i) = load(f_list(i).name, 'maxPeak');
end


meanDur_Append = [];
for k = 1:length(meanDur_all)
   meanDur_Append = [meanDur_Append; meanDur_all(k).meanDur];
end

maxPeak_Append = [];
for k = 1:length(maxPeak_all)
    maxPeak_Append = [maxPeak_Append; maxPeak_all(k).maxPeak];
end
        
cd(path);





% file_list = dir(path);
% 
% for i = 1:length(file_list);
%        if file_list(i).isdir && file_list(i).name(1) ~= '.'
%            cd([path filesep file_list(i).name]);
% %            matDir = dir('*.mat');
%            wvfFile = FindFiles('*.mat');
%            for j = 1:length(wvfFile) 
%                     if ~isempty(wvfFile)
%                         meanDur_all(j) = load(wvfFile{1}, 'meanDur');
%                         maxPeak_all(j) = load(wvfFile{1}, 'maxPeak');
%                     end
%            end
%        end
% end