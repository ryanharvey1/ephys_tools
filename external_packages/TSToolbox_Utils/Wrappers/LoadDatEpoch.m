function ep = LoadDatEpoch(varargin)
% 
% ep = LoadDatEpoch(epname)
% or
% ep = LoadDatEpoch(fbasename,epname)
% if you need to specify the path to merged folder

%Adrien Peyrache, 2012

if length(varargin)==1
 epname = varargin{1};
 fbasename = [];
elseif length(varargin)==2
 epname = varargin{2};
 fbasename = [varargin{1} filesep];
 
else
    error('bad number of arguments')
end
if exist([fbasename 'Analysis/BehavEpochs.mat'],'file');
    %warning off
    load([fbasename 'Analysis/BehavEpochs.mat'],epname);
    %warning on;
    if ~exist(epname)
        ep = intervalSet([],[]);
    else
        eval(['ep = ' epname ';'])
    end
else
    warning('LoadDatEpoch: no BehavEpochs file! Run CreateEpoch first')
    ep = intervalSet([],[]);
end
