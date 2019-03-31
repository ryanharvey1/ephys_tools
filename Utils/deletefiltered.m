function [counter]=deletefiltered(tfile,counter,event)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
counter = counter + 1;
[filepath, filename] = fileparts(tfile);
if exist([filepath filesep filename sprintf('S%d',event) '_spikeData.mat'],'file')
    delete([filepath filesep filename sprintf('S%d',event) '_spikeData.mat']);
elseif exist([filepath filesep filename '_spikeData.mat'],'file')
    delete([filepath filesep filename '_spikeData.mat']);
end
end

