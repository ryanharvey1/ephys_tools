function createChannelMapFile_Local(basepath)
%  create a channel map file based on the .xml file. Original script from
%  kilosort 1 wrapper. Modified to not use channels that user set as 'skip'
%  in the .xml file.
%
%  Modified by Eliezyer de Oliveira, 02/03/2020

if ~exist('basepath','var')
    basepath = cd;
end
d   = dir('*.xml');
par = LoadXml(fullfile(basepath,d(1).name));

xcoords = [];
ycoords = [];
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

for a= 1:ngroups %being super lazy and making this map with loops
    x = [];
    y = [];
    tchannels  = tgroups{a};
    for i =1:length(tchannels)
        x(i) = length(tchannels)-i;
        y(i) = -i*10;
        if mod(i,2)
            x(i) = -x(i);
        end
    end
    x = x+a*200;
    xcoords = cat(1,xcoords,x(:));
    ycoords = cat(1,ycoords,y(:));
end

Nchannels = length(xcoords);

kcoords = zeros(Nchannels,1);
for a= 1:ngroups
    kcoords(tgroups{a}+1) = a;
end

connected   = true(Nchannels, 1);
%getting channels labeled as a skip to set as disconnected
%with the code below I assume your anatomical group is the same as your
%spike group, so make sure it is if you want to skip channels correctly
aux_skip      = cell2mat({par.AnatGrps.Skip});
aux_ch        = cell2mat({par.AnatGrps.Channels});
disconnect_ch = aux_ch(logical(aux_skip))+1;

connected(disconnect_ch) = false;



chanMap     = 1:Nchannels;
chanMap0ind = chanMap - 1;
[~,I] =  sort(horzcat(tgroups{:}));
xcoords = xcoords(I);
ycoords  = ycoords(I);

save(fullfile(basepath,'chanMap.mat'), ...
    'chanMap','connected', 'xcoords', 'ycoords', 'kcoords', 'chanMap0ind')

%%
% % 
% Nchannels = 128;
% connected = true(Nchannels, 1);
% chanMap   = 1:Nchannels;
% chanMap0ind = chanMap - 1;
% 
% xcoords   = repmat([1 2 3 4]', 1, Nchannels/4);
% xcoords   = xcoords(:);
% ycoords   = repmat(1:Nchannels/4, 4, 1);
% ycoords   = ycoords(:);
% kcoords   = ones(Nchannels,1); % grouping of channels (i.e. tetrode groups)
% 
% save('C:\DATA\Spikes\Piroska\chanMap.mat', ...
%     'chanMap','connected', 'xcoords', 'ycoords', 'kcoords', 'chanMap0ind')
%%

% kcoords is used to forcefully restrict templates to channels in the same
% channel group. An option can be set in the master_file to allow a fraction 
% of all templates to span more channel groups, so that they can capture shared 
% noise across all channels. This option is

% ops.criterionNoiseChannels = 0.2; 

% if this number is less than 1, it will be treated as a fraction of the total number of clusters

% if this number is larger than 1, it will be treated as the "effective
% number" of channel groups at which to set the threshold. So if a template
% occupies more than this many channel groups, it will not be restricted to
% a single channel group. 