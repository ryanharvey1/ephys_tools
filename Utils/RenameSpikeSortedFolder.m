% RenameSpikeSortedFolder

% path to raw data
cd D:\Projects\PAE_PlaceCell\data

folders=dir('**/*SNAPSorterResults');
for i=1:size(folders,1)
    disp([folders(i).folder,filesep,folders(i).name])
    movefile([folders(i).folder,filesep,folders(i).name],...
        [folders(i).folder,filesep,'Sorted'])
end
disp('done')
