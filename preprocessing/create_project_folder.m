function create_project_folder(project_name)
% create_project_folder: sets up new project folder that conforms to 
% ephys_tools standards
%
% How to use:
% 1. cd to the location where you want the project to exist
% 2. type create_project_folder('project_name') with project_name replaced by
%       your project name into Command Window
% 3. Hit enter
% 4. Populate newly created folders following instructions below
%
%
% AnimalMetadata: should contain metadata for each animal. To create metadata
% files, run handleAnimalMetaData.m and follow Command Window instructions
%
% data: contains raw and spike sorted ephys data. Each animal should have
% their own folder labeled with their ID. Within each animal folder should
% be individual session folders from labeled with
% year-month-day_hour-minute-second (2019-08-15_19-25-40)
% For example
%   -HPCatn07
%       -2019-08-15_19-25-40
%           -CheetahLogFile
%           -CheetahLostADRecords
%           -CSC1
%           -CSC2
%           -CSC3
%           -CSC4
%           -TT1
%           -VT1
%
% ProcessedData: contains individual mat files for each session created 
% by postprocess.m. Naming scheme AnimalID_Syearmonthdayhourminutesecond 
% (HPCatn07_S20190815192540)
%
% Ryan Harvey 2019

mkdir(fullfile(pwd,project_name))
mkdir(fullfile(pwd,project_name,'AnimalMetadata'))
mkdir(fullfile(pwd,project_name,'data'))
mkdir(fullfile(pwd,project_name,'ProcessedData'))

end

