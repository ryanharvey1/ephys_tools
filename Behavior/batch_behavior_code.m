% Batch behavior_coding LB 09/June/2019
% runs a batch of directory names in random order 

path_to_files='G:\Maria\OF\Videos\lgOF_dlc';
files = dir([path_to_files,'\**\*.avi']);
% shuffVidIdx=randperm(size(files,1));

for i=1:size(files,1)
    path=fullfile(files(i).folder,filesep,files(i).name);
    behavior_coding('path',path,'behaviors',{'grooming','rearing','head_scan','Trial_start','Trial_end'}) 
    pause
    
end