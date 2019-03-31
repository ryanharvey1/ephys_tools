function PreprocessingMClust( input_,template,multible_sessions)
%PreprocessingMClust Summary of this function goes here
%   Detailed explanation goes here

addpath(genpath('F:\Users\BClarkLab\GoogleDrive\MatlabDir\MClust3.1UofL_c'));
addpath(genpath('d:\Users\BClarkLab\GoogleDrive\MatlabDir\MClust3.1UofL_c'));

if ~exist (input_,'dir')
    input_= uigetdir; 
end

Sep_ntt_by_event(input_)

output_ =[input_,'p'];

preprocessTTUofL(input_,output_,template,(1:8));

check_channels_kk(output_);


disp('write, save, and quit') % pause to write and okay saved batch file
for i=30:-1:1
    pause(1);
    disp([num2str(i),' Seconds Remaining'])
end

cd([output_ filesep 'TT']);
RunClustBatch(multible_sessions,'batch.txt');
% RunClustBatch2('batch.txt');

end

