% Mind_2017_Bayesian_decoding_testing
% https://nbviewer.jupyter.org/github/summer-mind/mind_2017/blob/master/Tutorials/SpikeDecoding/spike_decoding_matlab.ipynb#Bayesian-decoding
clear
% replace this with the location of your local GitHub root folder
SET_GitHub_root = '/Users/ryanharvey/GoogleDrive/MatlabDir';
% replace this with the location of your local data folder
% SET_data_fd = '/Users/ryanharvey/Downloads/R042-2013-08-18';
% warning is because this overloads the plot() function
addpath(genpath(cat(2,SET_GitHub_root,'/vandermeerlab/code-matlab/shared')));

%% Loading and inspecting the data
% Next, we load the spike and position data (stored in .t and .nvt files
% respectively, as elaborated on here, you may need to unzip the position data first):

data=load('/Users/ryanharvey/Downloads/HPCatn02_S20180706103807.mat');
% pos=[data.linear_track.nonlinearFrames(:,2)';data.linear_track.nonlinearFrames(:,3)'];
% ts=data.linear_track.nonlinearFrames(:,1)';
event=2;
pos=[data.frames(:,2)';data.frames(:,3)'];
ts=data.frames(:,1)';

pos=pos(:,ts>data.events(1,event) & ts<data.events(2,event));
ts=ts(:,ts>data.events(1,event) & ts<data.events(2,event));


remove=isnan(pos(1,:));
pos(:,remove)=[];
ts(:,remove)=[];


S=data.Spikes';
for i=1:length(S)
    S{i}=S{i}(S{i}>data.events(1,event) & S{i}<data.events(2,event));
end



%% Now, we can compute tuning curves for each of our cells.
% The TuningCurves() function does this in three steps:
% Compute an occupancy histogram (amount of time spent in each 2-D bin)
% Compute a spike histogram (count number of spikes fired in each 2-D bin)
% Compute firing rates by dividing the spike count by the occupancy
% The below code plots the output of the third step for our example cell:

nBinsx = round(100/3); nBinsy = round(100/3);
MinY = min(pos(2,:));
MaxY = max(pos(2,:));
MinX = min(pos(1,:));
MaxX = max(pos(1,:));
edges{1} = linspace(MinY, MaxY, nBinsy+1);
edges{2} = linspace(MinX, MaxX, nBinsx+1);

Omatrix = hist3([pos(2,:); pos(1,:)]','Edges',edges);

Omatrix(end,:) = [];
Omatrix(:,end) = [];
occ = Omatrix/30;
good_idx = find(Omatrix >= 1);


% bin spike data
for i=1:length(S)
    Smatrix = hist3([interp1(ts,pos(2,:),S{i}), interp1(ts,pos(1,:),S{i})],'Edges',edges);
    Smatrix(end,:) = [];
    Smatrix(:,end) = [];
    % divide binned spikes by occupancy to get rate maps and store data
    FilledRateMatrix = Smatrix./occ;
    FilledRateMatrix(isinf(FilledRateMatrix))=0;
    filtWidth = [5 5]; filtSigma = 1;
    imageFilter=fspecial('gaussian',filtWidth,filtSigma);
    SmoothRateMap(i,:,:) = nanconv(FilledRateMatrix,imageFilter, 'nanout');
    temp = SmoothRateMap(i,:,:);
    tc(i,:) = temp(good_idx);
    
end

%% You can inspect the tuning curves of all cells as follows:
ppf = 25; % plots per figure
for iC = 1:length(S)
    nFigure = ceil(iC/ppf); figure(nFigure);
    subtightplot(5,5,iC-(nFigure-1)*ppf);
    
    pcolor(sq(SmoothRateMap(iC,:,:))); shading flat; axis off;
    caxis([0 10]); % set color axis from 0 to 10 Hz
    
end
%% Obtaining the Q-matrix
cfg_Q = [];
cfg_Q.binsize = 0.25;
cfg_Q.tvec_edges = ts(1):cfg_Q.binsize:ts(end);
cfg_Q.tvec_centers = cfg_Q.tvec_edges(1:end-1)+cfg_Q.binsize/2;

for i=1:length(S)
    [Q.data(i,:),~] = histcounts(S{i},cfg_Q.tvec_edges);
end
% This ?Q-matrix? is of size [nCells x nTimeBins] and contains the spike
% count for each neuron in a given time bin:
% Q = MakeQfromS(cfg_Q,DEC_S);
figure;
imagesc(Q.data);colormap( flipud(gray))
%% Bayesian decoding
% As noted in the introduction above, given that we have neurons whose
% activity seems to encode some stimulus variable (location in this case,
% as evident from the tuning curves above), we can attempt to decode that
% variable based on the neurons' time-varying activity.
% A popular approach to doing this is ?one-step Bayesian decoding?,
% illustrated in this figure (from van der Meer et al. 2010):


nBins = length(tc);
occUniform = repmat(1/nBins,[nBins 1]); % prior over locations, P(x) in section above

nActiveNeurons = sum(Q.data > 0);

len = size(Q.data,2);
p = nan(size(Q.data,2),nBins); % this variable will store the posterior
for iB = 1:nBins % loop over space bins (x_i)
    % these 3 lines implement the actual decoding computation
    tempProd = nansum(log(repmat(tc(:,iB),1,len).^Q.data));
    % compare to the equations above!
    tempSum = exp(-cfg_Q.binsize*nansum(tc(:,iB)',2));
    p(:,iB) = exp(tempProd)*tempSum*occUniform(iB);
end

p = p./repmat(sum(p,2),1,nBins); % renormalize to 1 total probability
p(nActiveNeurons < 1,:) = 0; % ignore time bins with no activity
%% Now we want to display the results. Before we do so, we should convert
% the rat's actual position into our binned form, so that we can compare it
% to the decoded estimate:

xBinned = interp1(ts,rescale(pos(1,:),1,nBinsx),cfg_Q.tvec_edges);
yBinned = interp1(ts,rescale(pos(2,:),1,nBinsy),cfg_Q.tvec_edges);

%%
dec_err = nan(length(Q.data),1); % variable to keep track of decoding error
h = figure; set(h,'Position',[100 100 640 480]); % fix position of figure to help encoding of movie for export later

for iT = 1:length(Q.data) % loop over time bins
    
    toPlot = nan(nBinsy,nBinsx); % initialize empty 2-D grid
    toPlot(good_idx) = p(iT,:); % assign decoding probabilities to "good" bins
    
    toPlot = flipud(rot90(toPlot));
    
    cla; pcolor(toPlot); axis xy; hold on; caxis([0 0.5]); shading flat; axis off; % plot decoding output
    hold on; plot(yBinned(iT),xBinned(iT),'.r','MarkerSize',30); % overlay actual position
    
    % get x and y coordinates of location with largest decoding probability ("maximum a posteriori" or MAP)
    [~,idx] = max(toPlot(:)); [x_map,y_map] = ind2sub(size(toPlot),idx);
    plot(y_map,x_map,'g*','MarkerSize',5); % plot MAP location
    decode_x(iT)=x_map;
    decode_y(iT)=y_map;
    
    % compute error: distance between actual position and MAP position
    if nActiveNeurons(iT) > 0, dec_err(iT) = sqrt((yBinned(iT)-y_map).^2+(xBinned(iT)-x_map).^2); end
    
    h = title(sprintf('t %.2f, nCells %d, dist %.2f',cfg_Q.tvec_edges(iT),nActiveNeurons(iT),dec_err(iT)));
    if nActiveNeurons(iT) == 0, set(h,'Color',[1 0 0]); else, set(h,'Color',[0 0 0]); end
    
%         pause(0.1); drawnow;
    % f(iT) = getframe(gcf); % uncomment this if exporting to video file
    
end

figure;
plot(xBinned,yBinned,'.k')
figure;
plot(decode_x,decode_y,'.r')

figure;
scatter(xBinned(1:end-1),decode_x)
figure;
scatter(yBinned(1:end-1),decode_y)
%%
% dec_err_tsd = tsd(cfg_Q.tvec_edges,dec_err); % create timestamped data struct for decoding error computed above
% 
% cfg = []; cfg.y_edges = nBinsy; cfg.x_edges = nBinsx;
% Enc.x=pos(1,:);
% Enc.y=pos(2,:);
% space_err = TSDbySpace(cfg,Enc,dec_err_tsd);
% 
% figure;
% pcolor(space_err); shading flat; axis off; colorbar; caxis([0 10]);




%%