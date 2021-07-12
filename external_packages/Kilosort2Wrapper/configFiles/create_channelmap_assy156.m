
function create_channelmap_assy156(basepath,varargin)
%  create a channel map file based on the .xml file. Original script from
%  kilosort 1 wrapper. Modified to not use channels that user set as 'skip'
%  in the .xml file. Modified further to accomodate probe geometry. 
%
%  Default map is for ASSY-156-E1 64 channel Cambridge Neurotech Probe
%
%  Modified by Eliezyer de Oliveira, 02/03/2020
%  Modified by Laura Berkowitz, 11/29/2020 

p = inputParser; 
p.addParameter('xcoords',[0,70,5,65,...% shank A
    10,60,15,55,...
    20,50,25,45,...
    30,40,35,35,...
    0+250,70+250,5+250,65+250,... % shank B
    10+250,60+250,15+250,55+250,...
    20+250,50+250,25+250,45+250,...
    30+250,40+250,35+250,35+250,...
    0+500,70+500,5+500,65+500,... % shank C
    10+500,60+500,15+500,55+500,...
    20+500,50+500,25+500,45+500,...
    30+500,40+500,35+500,35+500,...
    0+750,70+750,5+750,65+750,... % shank D
    10+750,60+750,15+750,55+750,...
    20+750,50+750,25+750,45+750,...
    30+750,40+750,35+750,35+750]');
p.addParameter('ycoords', [0,-25,-40,-65,...% shank A
    -80, -105, -120, -145,...
    -160,-185,-200,-225,...
    -240,-265,-280,-305,...
    0,-25,-40,-65,... % shank B
    -80, -105, -120, -145,...
    -160,-185,-200,-225,...
    -240,-265,-280,-305,...
    0,-25,-40,-65,... % shank C
    -80, -105, -120, -145,...
    -160,-185,-200,-225,...
    -240,-265,-280,-305,...
    0,-25,-40,-65,... % shank D
    -80, -105, -120, -145,...
    -160,-185,-200,-225,...
    -240,-265,-280,-305]');
p.addParameter('Nshanks',4)
p.addParameter('fs',30000)
p.parse(varargin{:});
if ~exist('basepath','var')
    basepath = cd;
end

d   = dir('*.xml');
par = LoadXml(fullfile(basepath,d(1).name));

xcoords = p.Results.xcoords;
ycoords = p.Results.ycoords;
if ~isfield(par,'nElecGps')
    warning('No Electrode/Spike Groups found in xml.  Using Anatomy Groups instead.')
    tgroups = par.ElecGp;
    ngroups = length(tgroups);
else
    t = par.AnatGrps;
    ngroups = length(par.AnatGrps);
    for g = 1:ngroups
        tgroups{g} = par.AnatGrps(g).Channels;
    end
end

% Unpack to save
Nchannels = length(xcoords);
Nshanks = p.Results.Nshanks;
fs = p.Results.fs;

kcoords = (reshape(repmat(1:Nshanks, Nchannels/Nshanks, 1), Nchannels, 1));

connected   = true(Nchannels, 1); 
%getting channels labeled as a skip to set as disconnected
%with the code below I assume your anatomical group is the same as your
%spike group, so make sure it is if you want to skip channels correctly
aux_skip      = cell2mat({par.AnatGrps.Skip});
aux_ch        = cell2mat({par.AnatGrps.Channels});
disconnect_ch = aux_ch(logical(aux_skip))+1;

chanMap   = 1:Nchannels;
chanMap0ind = chanMap - 1;

connected(disconnect_ch) = false;

save(fullfile(basepath,'chanMap.mat'), ...
    'chanMap','connected', 'xcoords', 'ycoords', 'kcoords', 'chanMap0ind', 'fs')
end