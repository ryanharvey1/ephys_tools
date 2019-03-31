
%% matlab search path
home = pwd;

addpath(home);
addpath([home '\Neuralynx']);
addpath([home '\Utils']);
addpath([home '\Classes']);
addpath([home '\IO']);
addpath([home '\MathStats']);
addpath([home '\Analysis']);

%% other settings
more;                % set display page mode ON
cd(home);            % set current working directory

fprintf(2, 'NSMA\\Startup.m: NSMA-2 path loaded!\n');
