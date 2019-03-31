function [returnVar,msg] = RemoveDCfromDat_AllDat(fbasename)

% USAGE:
%     RemoveDCfromDat_AllDat(fbasename)
% 
% removes DC offset in all channels from all dat files starting with
% 'fbasename'. Assumes that there is at least one xml file associated to
% one the dat file (must be created before with Neuroscope and/or NDmanger)
%     
% Adrien Peyrache, 2012

d = dir([fbasename '*.dat']);
ii = 1;
fname = d(ii).name;
fname = fname(1:end-4);
while ~exist([fname '.xml'],'file')
    fname = d(ii).name;
    fname = fname(1:end-4);
    ii=ii+1;
end
if exist([fname '.xml'],'file')
    syst = LoadXml(fname);
else
    error('No XML file!')
end

for ii=1:length(d)
    [returnVar,msg] = RemoveDCfromDat(d(ii).name,syst.nChannels);
    if ~returnVar
        warning(msg)
    else
        printf('...done\n');
    end
end