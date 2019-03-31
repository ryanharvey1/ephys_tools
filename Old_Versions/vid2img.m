%% This script is to be used when converting video files to image sequences (for use in imageJ). 
% Laura B & Ryan H  June 2016

disp('locate video folder');
workingDir= uigetdir;                                           % select folder containing videos
MClustPath = '/Users/RyanHarvey/Dropbox/MATLAB/Cell Analysis/'; % adds function 'findfiles.m' to path
    if exist('FindFiles.m','dir')<0;                            % Opens up file window for you to locate 'FindFiles'
            disp('locate FindFiles'); 
            MClustPath = uigetdir;    
    end
addpath(genpath(MClustPath));                                   % you can do this manually before running script
videofile = FindFiles('*.MOV','StartingDirectory', workingDir); % locates how many .mov files are in selected folder
counter = 1;
for i = 1:length(videofile);                                    % starts loop for number of .mov files available in folder      
      name=dir(videofile{i});                                   % Locates video file name
        for k=1:length(name);
            videofiles=name(k).name;                            % extracts video file name from structure
        end
    mwmVideo = VideoReader(videofiles);                         % Create VideoReader - create a VideoReader to use for reading frames from
    ii=1;                                                       % Create the image sequence
    [filepath, filename] = fileparts(videofile{i});
    mkdir([filepath filesep filename 'images']);
    foldername = [filepath filesep filename 'images'];
        while hasFrame(mwmVideo)
            img = readFrame(mwmVideo);
            filename=[sprintf('%03d',ii) '.jpg'];
            fullname=fullfile(foldername,filename);
            imwrite(img,fullname)                               % Write out to a JPEG file (img1.jpg, img2.jpg, etc.)
            ii = ii+1;
        end
    fprintf('Just finished iteration #%d\n', counter);
    counter = counter + 1;
end
% Remember to clear all before running for a second time ***

