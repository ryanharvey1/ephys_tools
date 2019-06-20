% getCueCoords finds boundary coordinates of maze objects.
% LB May 2019
%   ASSUMPTIONS: 
%    Videos must be in format compatible with Video Reader (see matlab
%   documentation for compatibilities given your computers OS). Recommend .AVI
%   as its compatible for most OS.
%
%    Assumes videos and images are stored in subdirectories. Saves first frame used to
%   grab XY in same location as data video or grabs image for obtaining
%   coordinates (should be same name as subdirectory).
%   
%   Subdirectories names become subID in saved table. 

%   INPUT: 
%       vidDir --> path to video directory in string format e.g.
%       'G:\Maria\OF\Videos\lgOF_dlc'
%   OUTPUT: 
%       cueCoords --> table containing x y coordinates for the object boundary.  
%
% 

function cueCoords=getCueCoords(vidDir)
vidfile=dir(vidDir);
vidfile(1:2,:)=[]; %get rid of . & .. entries
cueCoords=table; %initialize data table.

for file=1:length(vidfile) %loop through folders containing subject videos
    %Pull up video
    videoPathAndFileName = strcat(vidfile(file).folder,filesep,vidfile(file).name,filesep,vidfile(file).name,'.avi');
    videoObj = VideoReader(videoPathAndFileName);
    
    %Save first frame as image (jpeg)
    if isempty(dir('*.jpeg'))
        imwrite(readFrame(videoObj),[vidfile(file).folder,filesep,vidfile(file).name,filesep,vidfile(file).name,'.jpeg'])
        cd(strcat(vidfile(file).folder,filesep,vidfile(file).name))% go to folder containing video
        imshow(readFrame(videoObj))
        hold on
    else
        cd(strcat(vidfile(file).folder,filesep,vidfile(file).name))% go to folder containing image
        imshow(imread([vidfile(file).name, '.jpeg'])) %display the first frame
        hold on;
    end
    
    answer = questdlg('Would you like to track a cue?', ...
        'Cue Check', ...
        'Yes','No','No');
    % Handle response
    switch answer
        case 'Yes'
            disp([answer ' let''s get some coordinates!'])
            leave=0;
        case 'No'
            disp([answer ' on to the next one!'])
            leave=1;
    end
    
    if leave==1
         cueCoords.subID{file}=vidfile(file).name; %save the id
        cueCoords.coords{file}=NaN; 
        continue
    end
    
    %set iteration
    i=1;
    % let the user click around the boundary (code from Ryan H)
    while true
        [X,Y]=ginput(1);
        if isempty(X)
            break
        end
        corners(i,:)=[X,Y];
        plot(corners(:,1),corners(:,2),'r',X,Y,'*r')
        i=i+1;
    end
    
    cueCoords.subID{file}=vidfile(file).name; %save the id
    cueCoords.coords{file}=corners; 
    
    clear corners
end
end