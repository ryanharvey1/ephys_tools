% testingFindLoops
clear;clc;close all
addpath('/Users/ryanharvey/GoogleDrive/MatlabDir/BClarkToolbox/Analysis/MWM')
cd ('/Volumes/Samsung USB/RyanData')
files=dir;
for i=20:length(files)
%     disp([num2str(i),' of ',num2str(length(files))])
    XY=csvread(files(i).name,5,2);
    figure(1);plot(XY(:,1),XY(:,2));hold on;pause(.0000001)
    FinalLoop=FindLoops(XY(:,1),XY(:,2));
    size(FinalLoop)
    try
        h=streamline({FinalLoop.loops});
        set(h,'Color','red');
        pause(5);
    catch
    end
    close all
    clear FinalLoop
end