% Mind_2017_Bayesian_decoding
% https://nbviewer.jupyter.org/github/summer-mind/mind_2017/blob/master/Tutorials/SpikeDecoding/spike_decoding_matlab.ipynb#Bayesian-decoding
clear
 % replace this with the location of your local GitHub root folder
SET_GitHub_root = '/Users/ryanharvey/GoogleDrive/MatlabDir';
% replace this with the location of your local data folder
SET_data_fd = '/Users/ryanharvey/Downloads/R042-2013-08-18'; 
% warning is because this overloads the plot() function
addpath(genpath(cat(2,SET_GitHub_root,'/vandermeerlab/code-matlab/shared'))); 

%% Loading and inspecting the data
% Next, we load the spike and position data (stored in .t and .nvt files
% respectively, as elaborated on here, you may need to unzip the position data first):
cd(SET_data_fd);
 
please = []; please.load_questionable_cells = 1;
S = LoadSpikes(please); % `please` variable overrides default LoadSpikes() options
 
pos = LoadPos([]); % empty input [] causes LoadPos() to use default options
%% Before looking at the data, we will first exclude the pre- and post-run 
% segments of the data:

LoadExpKeys; % annotation file containing some basic information about this data session
S = restrict(S,ExpKeys.TimeOnTrack,ExpKeys.TimeOffTrack); % restrict to on-track data only
pos = restrict(pos,ExpKeys.TimeOnTrack,ExpKeys.TimeOffTrack);
%% Now we can plot the position data (location of the rat as it runs the 
% maze, in gray), and the position of the rat when an example neuron fired 
% a spike (in red):
plot(getd(pos,'y'),getd(pos,'x'),'.','Color',[0.5 0.5 0.5],'MarkerSize',1); % note getd() gets data corresponding to a specific label (x and y here)
axis off; hold on;

iC = 7; % select cell 7 (out of 107 total)
spk_x = interp1(pos.tvec,getd(pos,'x'),S.t{iC},'linear');
spk_y = interp1(pos.tvec,getd(pos,'y'),S.t{iC},'linear');
 
h = plot(spk_y,spk_x,'.r'); axis tight;
%% Estimating tuning curves
% This figure is a useful visualization of the raw data, but it is not a 
% tuning curve, which captures the relationship between a variable of 
% interest (e.g. position) to firing rate. A set of tuning curves can be 
% thought of as an encoding model that specifies how the position variable 
% is encoded by our population of place cells. As a first step to 
% estmiating this model, we restrict our data to those times the rat is 
% actually running on the track:
LoadMetadata; % loads experimenter-annotated file associated with each data session

% ENCoding variables: used for estimating tuning curves (encoding model)
ENC_S = restrict(S,metadata.taskvars.trial_iv); % trial_iv contains the start and end times of trials
ENC_pos = restrict(pos,metadata.taskvars.trial_iv);
 
% check for empties and remove
keep = ~cellfun(@isempty,ENC_S.t);
ENC_S = SelectTS([],ENC_S,keep);

% also set up DECoding variables for use later
DEC_S = SelectTS([],S,keep);
%% Now, we can compute tuning curves for each of our cells. 
% The TuningCurves() function does this in three steps:
% Compute an occupancy histogram (amount of time spent in each 2-D bin)
% Compute a spike histogram (count number of spikes fired in each 2-D bin)
% Compute firing rates by dividing the spike count by the occupancy
% The below code plots the output of the third step for our example cell:
cfg_tc = [];
cfg_tc.minOcc = 1; % minimum occupancy (in seconds) for bin to be included
cfg_tc.binEdges{1} = 80:10:660; cfg_tc.binEdges{2} = 0:10:520; % bin edges (in camera pixels)
tc = TuningCurves(cfg_tc,ENC_S,ENC_pos);

pcolor(sq(tc.tc2D(iC,:,:))); % sq() squeezes 3-D matrix down to 2-D for plotting
shading flat; axis off; colorbar;
%%
cfg_tc.smoothingKernel = gausskernel([4 4],2); % Gaussian kernel of 4x4 pixels, SD of 2 pixels (note this should sum to 1)
tc = TuningCurves(cfg_tc,ENC_S,ENC_pos);
figure
pcolor(sq(tc.tc2D(iC,:,:))); shading flat; axis off; colorbar
%% You can inspect the tuning curves of all cells as follows:
ppf = 25; % plots per figure
for iC = 1:length(ENC_S.t)
    nFigure = ceil(iC/ppf); figure(nFigure);
    subtightplot(5,5,iC-(nFigure-1)*ppf);
    
    pcolor(sq(tc.tc2D(iC,:,:))); shading flat; axis off;
    caxis([0 10]); % set color axis from 0 to 10 Hz
 
end
%% Obtaining the Q-matrix
cfg_Q = [];
cfg_Q.binsize = 0.25;
cfg_Q.tvec_edges = metadata.taskvars.trial_iv.tstart(1):cfg_Q.binsize:metadata.taskvars.trial_iv.tend(end);
cfg_Q.tvec_centers = cfg_Q.tvec_edges(1:end-1)+cfg_Q.binsize/2;

% This ?Q-matrix? is of size [nCells x nTimeBins] and contains the spike 
% count for each neuron in a given time bin:
Q = MakeQfromS(cfg_Q,DEC_S);
figure;
imagesc(Q.data);colormap( flipud(gray))
%% Bayesian decoding
% As noted in the introduction above, given that we have neurons whose 
% activity seems to encode some stimulus variable (location in this case, 
% as evident from the tuning curves above), we can attempt to decode that 
% variable based on the neurons' time-varying activity.
% A popular approach to doing this is ?one-step Bayesian decoding?, 
% illustrated in this figure (from van der Meer et al. 2010):

Q = restrict(Q,metadata.taskvars.trial_iv); % for speed, only decode trials of running on track

nBins = length(tc.usr.good_idx);
occUniform = repmat(1/nBins,[nBins 1]); % prior over locations, P(x) in section above

nActiveNeurons = sum(Q.data > 0);
 
len = length(Q.tvec);
p = nan(length(Q.tvec),nBins); % this variable will store the posterior
for iB = 1:nBins % loop over space bins (x_i)
    % these 3 lines implement the actual decoding computation
    tempProd = nansum(log(repmat(tc.tc(:,iB),1,len).^Q.data));
    % compare to the equations above!
    tempSum = exp(-cfg_Q.binsize*nansum(tc.tc(:,iB)',2)); 
    p(:,iB) = exp(tempProd)*tempSum*occUniform(iB);
end
 
p = p./repmat(sum(p,2),1,nBins); % renormalize to 1 total probability
p(nActiveNeurons < 1,:) = 0; % ignore time bins with no activity
%% Now we want to display the results. Before we do so, we should convert 
% the rat's actual position into our binned form, so that we can compare it
% to the decoded estimate:
xBinned = interp1(ENC_pos.tvec,tc.usr.pos_idx(:,1),Q.tvec);
yBinned = interp1(ENC_pos.tvec,tc.usr.pos_idx(:,2),Q.tvec);
%%
dec_err = nan(length(Q.tvec),1); % variable to keep track of decoding error
h = figure; set(h,'Position',[100 100 640 480]); % fix position of figure to help encoding of movie for export later

for iT = 1:length(Q.tvec) % loop over time bins
    
    toPlot = nan(tc.usr.nBins{1},tc.usr.nBins{2}); % initialize empty 2-D grid
    toPlot(tc.usr.good_idx) = p(iT,:); % assign decoding probabilities to "good" bins
 
    cla; pcolor(toPlot); axis xy; hold on; caxis([0 0.5]); shading flat; axis off; % plot decoding output
    hold on; plot(yBinned(iT),xBinned(iT),'ow','MarkerSize',15); % overlay actual position
 
    % get x and y coordinates of location with largest decoding probability ("maximum a posteriori" or MAP)
    [~,idx] = max(toPlot(:)); [x_map,y_map] = ind2sub(size(toPlot),idx);
    % plot(y_map,x_map,'g*','MarkerSize',5); % plot MAP location
 
    % compute error: distance between actual position and MAP position
    if nActiveNeurons(iT) > 0, dec_err(iT) = sqrt((yBinned(iT)-y_map).^2+(xBinned(iT)-x_map).^2); end
 
    h = title(sprintf('t %.2f, nCells %d, dist %.2f',Q.tvec(iT),nActiveNeurons(iT),dec_err(iT))); 
    if nActiveNeurons(iT) == 0, set(h,'Color',[1 0 0]); else, set(h,'Color',[0 0 0]); end
    
    pause(0.1); drawnow;
    % f(iT) = getframe(gcf); % uncomment this if exporting to video file
   
end
%%
dec_err_tsd = tsd(Q.tvec,dec_err); % create timestamped data struct for decoding error computed above

cfg = []; cfg.y_edges = cfg_tc.binEdges{2}; cfg.x_edges = cfg_tc.binEdges{1};
space_err = TSDbySpace(cfg,ENC_pos,dec_err_tsd);
 
figure;
pcolor(space_err); shading flat; axis off; colorbar; caxis([0 10]);
%%