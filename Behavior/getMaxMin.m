% getMaxMin finds max/min values of maze apparatus from first frame of video files.
% LB May 2019
%   ASSUMPTIONS: 
%    Videos must be in format compatible with Video Reader (see matlab
%   documentation for compatibilities given your computers OS). Recommend .AVI
%   as its compatible for most OS.
%
%    Assumes videos are stored in subdirectories. Saves first frame used to
%   grab XY in same location as data video.
%   
%   Subdirectories names become subID in saved table. 
%   INPUT: 
%       vidDir --> path to video directory in string format e.g.
%       'G:\Maria\OF\Videos\lgOF_dlc'
%   OUTPUT: 
%       maxmin --> table containing max min values for each video file. 
%
% saves

function maxmin=getMaxMin(vidDir)
vidfile=dir(vidDir);
vidfile(1:2,:)=[]; %get rid of . & .. entries
maxmin=table; %initialize data table.

for file=1:length(vidfile) %loop through folders containing subject videos
    %Pull up video
    videoPathAndFileName = strcat(vidfile(file).folder,filesep,vidfile(file).name,filesep,vidfile(file).name,'.avi');
    videoObj = VideoReader(videoPathAndFileName);
    
    %Save first frame as image (jpeg)
    imwrite(readFrame(videoObj),[vidfile(file).folder,filesep,vidfile(file).name,filesep,vidfile(file).name,'.jpeg'])
    
    % go to folder containing video & image
    cd(strcat(vidfile(file).folder,filesep,vidfile(file).name))
    imshow(readFrame(videoObj)) %display the first frame
    hold on
    i=1;
    
    % let the user click around the coordinates
    while true
        [X,Y]=ginput(1);
        if isempty(X)
            break
        end
        corners(i,:)=[X,Y];
        plot(corners(:,1),corners(:,2),'r',X,Y,'*r')
        i=i+1;
    end
    
    maxmin.subID{file}=vidfile(file).name;
    maxmin.xmax{file}=max(corners(:,1)); maxmin.xmin{file}=min(corners(:,1));
    maxmin.ymax{file}=max(corners(:,2)); maxmin.ymin{file}=min(corners(:,2));
    
    clear corners
end
end