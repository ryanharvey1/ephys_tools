function [corrected_p,corrected_d,fit]=FMLM_wrapper(data,cell,session,varargin)
% FMLM_wrapper: wrapper to use the factorial Maximum Likelihood model with 
% the ephys_tools data structure 'data'

% parse inputs
p = inputParser;
p.addParameter('fig',0);
p.addParameter('dir_bin_width',6);
p.addParameter('pos_smooth',3);
p.addParameter('dir_smooth',1);
p.addParameter('max_iter',2000);
p.addParameter('accuracy',0.0001);
p.addParameter('tol',0.1);
p.parse(varargin{:});

fig = p.Results.fig;
dir_bin_width = p.Results.dir_bin_width;
pos_smooth = p.Results.pos_smooth;
dir_smooth = p.Results.dir_smooth;
max_iter = p.Results.max_iter;
accuracy = p.Results.accuracy;
tol = p.Results.tol;

% get spike binary matrix
[data_video_spk,occframes]=createframes_w_spikebinary(data,session,cell);
spikeframes=data_video_spk(data_video_spk(:,6)==1,:);

% find xy edges
nBinsx = round(data.maze_size_cm(session)/3); 
nBinsy = round(data.maze_size_cm(session)/3);
MinY = min(occframes(:,3));
MaxY = max(occframes(:,3));
MinX = min(occframes(:,2));
MaxX = max(occframes(:,2));
edges{1} = linspace(MinY, MaxY, nBinsy+1);
edges{2} = linspace(MinX, MaxX, nBinsx+1);

% loop through directions to create time and spike matrix
for a=0:dir_bin_width:360-6
    occ=occframes(occframes(:,4)>a & occframes(:,4)<a+6,:);
    occ_bin = hist3([occ(:,3) occ(:,2)],'Edges',edges);
    occ_bin(end,:) = [];
    occ_bin(:,end) = [];
    times(a/dir_bin_width+1,:,:) = occ_bin/data.samplerate;
    
    spk=spikeframes(spikeframes(:,4)>a & spikeframes(:,4)<a+6,:);
    spk_bin = hist3([spk(:,3) spk(:,2)],'Edges',edges);
    spk_bin(end,:) = [];
    spk_bin(:,end) = [];
    spikes(a/dir_bin_width+1,:,:) = spk_bin;
end

if fig
    figu=figure('Name',[data.rat,'  ',data.sessionID,'  ',...
        data.spikesID.TetrodeNum{cell},' Cell: ',...
        num2str(data.spikesID.CellNum(cell))],'NumberTitle','off');
    figu.Color=[1 1 1];
end

% get corrected position and directional tuning curves
[corrected_p, corrected_d, uncorrected_p, uncorrected_d,fit] = pxd(spikes,...
    times, pos_smooth, dir_smooth, max_iter, accuracy, tol, fig);
end