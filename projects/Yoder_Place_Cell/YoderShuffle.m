% YoderShuffle
clear; close all; clc
% data=importdata('/Users/RyanHarvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/GW- Field_StatsPlace Cells_Tilted_Mice_1_1.xlsx');
load('/Users/RyanHarvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/GW_Field_StatsPlaceCells_Tilted_Mice.mat')
placecells=data.textdata.Field_StatsPlaceCells0x2DTilted(2:end,1);

for i=1:2
    S1_S2(:,:,i)=placecells(i:5:end,:);
end

rSHUFF=[];
for i=1:length(S1_S2(:,:,2))
    [path,name,ext]=fileparts(S1_S2{i,1,2});
    load(strcat(path,filesep,name,'_Data.mat'),'field','SaveField');
    disp(['SHUFFLING:',strcat(path,filesep,name,'_Data.mat')])
    if exist('field','var') && exist('SaveField','var')
        for x = 1:1000
            field(randperm(numel(field)))=field;
            [D,M]=Displacement(SaveField,field);
            rSHUFF=[rSHUFF;M];
        end
    end
    clear SaveField field
end

NintyFith_Percentile = prctile(rSHUFF,95)
