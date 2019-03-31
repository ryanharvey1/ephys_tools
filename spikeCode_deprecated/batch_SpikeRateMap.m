% oneCell_SpikesRateMap creates summary of data in form of spikes on path and rate
% map for place cell and grid cell analysis
%
% Input: 
%       - raw data output from Labview
%
% Output:
%       - rate map stats (sparsity, information content)
%       - SpatialAutoCorr = spatial autocorrelation
%       - SmoothRateMap = smoothed rate map
%       - SessStability = within session stability
%       - GridScore
% 
% Created by Ben C Jan 2014

% identify path to data files
path = 'C:\Users\ben.clark\Dropbox\EC_analysis\__RawData_&_Summary\__control_results_files_only';
ReadData = FindFiles('*.txt', 'StartingDirectory', path);

for i = 1:length(ReadData);
    % video frame rate in Hz
    sampleRate = 60;

    % maze parameters for non-detects
    xmin = 0;
    xmax = 255;
    ymin = 0;
    ymax = 255;

    % min and max LED distance parameters
    minLED = 4;
    maxLED = 39;

    % load xy data
    data = importdata(ReadData{i});

    % extract coords, spikes, angle, and direction from ReadData output
    rx = data.data(:,2); % red x-coords
    ry = data.data(:,3); % red y-coords
    gx = data.data(:,4); % green x-coords
    gy = data.data(:,5); % green y-coords
    spks = data.data(:,6); % spikes
    angle = data.data(:,10); % angle in radians
    distance = data.data(:,11); % distance in pixels between LEDs
    datai = [rx,ry,gx,gy,spks,angle,distance]; % create array with all variables

    % find red LED non-detects
    dataFiltx = find(datai(:,1) > xmin & datai(:,1) < xmax & datai(:,2) < ymax & datai(:,1) > ymin);
    rFILT = datai(dataFiltx,:);

    % find green LED non-detects
    dataFiltxy = find(rFILT(:,3) > xmin & rFILT(:,3) < xmax & rFILT(:,4) < ymax & rFILT(:,4) > ymin);
    rgFILT = rFILT(dataFiltxy,:);

    % find Min and Max LED distance
    dataFiltxLED = find(rgFILT(:,7) > minLED & rgFILT(:,7) < maxLED);
    rgmmFILT = rgFILT(dataFiltxLED,:);

    % extract x y coords and smooth 
    rxs = runline(rgmmFILT(:,1),5,1); % smooth with 5pt window and 1pt step (from Chronux toolbox)
    rys = runline(rgmmFILT(:,2),5,1);

    % velocity of rat from smoothed xy data
    vel_x = diff(rxs); % vel units are pixels/frame
    vel_y = diff(rys);
    vel_abs = sqrt(vel_x.^2 + vel_y.^2); % scalar length of velocity vector = "scalar velocity" in pixels/frame
    vel_cmPerSec = vel_abs * 2.1 * sampleRate; % ~2.1 cm/pixel for Taube lab according to Shawn W  

    % find xy coords and spikes when rat is above velocity threshold
    iVel = find(vel_cmPerSec >= 0); % 5cm/sec is used by Stensola et al 2012
    spksf = rgmmFILT(:,5);
    dataVEL = [rxs rys spksf];
    dataVEL = dataVEL(iVel,:);

    % extract spike locations from smoothed and velocity filtered data
    spksi = find(dataVEL(:,3) == 1); 
    spks_VEL = dataVEL(spksi,:);

    % create array for binning data
    nBins = 26; % adjusted so the bins are ~5x5cm as in Koenig et al 2011
    MinY = min(dataVEL(:,2));
    MaxY = max(dataVEL(:,2));
    MinX = min(dataVEL(:,1));
    MaxX = max(dataVEL(:,1));
    edges{1} = linspace(MinY, MaxY, nBins+1);
    edges{2} = linspace(MinX, MaxX, nBins+1);

    % bin occupancy data
    occMatrix = [dataVEL(:,2),dataVEL(:,1)];
    Omatrix = hist3(occMatrix,'Edges',edges);
    Omatrix(1,:) = [];
    Omatrix(:,end) = [];
    occ = Omatrix/sampleRate;

    % bin spike data
    spikeMatrix = [spks_VEL(:,2), spks_VEL(:,1)];
    Smatrix = hist3(spikeMatrix,'Edges',edges);
    Smatrix(1,:) = [];
    Smatrix(:,end) = [];

    % divide binned spikes by occupancy to get rate maps and store data
    % (can use occ+eps instead of removing 0's)
    BFiringRateMatrix = Smatrix./occ;
    FilledRateMatrix = BFiringRateMatrix;
    FilledRateMatrix(isnan(FilledRateMatrix)) = 0;
    FilledRateMatrix(isinf(FilledRateMatrix)) = 0;

    % smooth raw rate histogram using a guassian (from NSMA toolbox)
    [SmoothRateMap] = smooth(FilledRateMatrix,1,5,1,5); % smooth with 5x5pt window and 1pt step

    % calculate rate map statistics
    [SessStability,FirstHalfSess,SecHalfSession] = SessionStability(rgmmFILT,nBins,sampleRate); % within-session stability
    rY = reshape(SmoothRateMap,nBins*nBins,1); % reshape data into column
    occRSHP = reshape(occ,nBins*nBins,1); % reshape data into column
    occSUM = sum(occRSHP); % summed occupancy
    pX = occRSHP./occSUM; % normalized occupancy
    [InformationContent] = InformationPerSpike(rY,pX); % information content from NSMA toolbox
    [Sparsity] = Sparsity(rY'); % sparsity from NSMA toolbox

    % create spatial autocorr and calculate gridscore 
    [SpatialAutoCorr] = SpatialAutoCorr(SmoothRateMap,nBins);
    [maxGS1, maxGS2, GS1, GS2, rGridScore, centerCUT] = GridScore(SpatialAutoCorr,nBins);

    % plot smoothed spike on path plot using red LED
    figure (1), plot(dataVEL(:,1),dataVEL(:,2),'LineWidth',1,'color',[0,0,0]+0.8); 
    axis square tight
    hold on
    scatter(spks_VEL(:,1),spks_VEL(:,2),45,'filled','k');
    box off
    title('Spike (black dots) on Path (gray)');
    print(figure (1), '-djpeg', filename);
    close(figure (1))
    
    % plot filled and smoothed rate map 
    figure(2), h = pcolor(SmoothRateMap);
    axis square tight
    hold on
    colorbar;
    box off
    set(h, 'EdgeColor', 'none');
    title('Smoothed Rate Map Entire Session');
    
    % plot spatial autocorr and plot with center cut out
    figure(3), h = pcolor(SpatialAutoCorr);
    axis square tight
    hold on
    colorbar;
    box off
    set(h, 'EdgeColor', 'none');
    title('Spatial Autocorrelation');
    
    % plot spatial auto corr with center and outside of ring removed
    figure(4), h = pcolor(centerCUT);
    axis square tight
    hold on
    colorbar;
    box off
    set(h, 'EdgeColor', 'none');

    % data and images to retain
    keep('maxGS1', 'maxGS2', 'GS1', 'GS2', 'rGridScore', 'SpatialAutoCorr', 'SmoothRateMap', 'SessStability', 'FirstHalfSess', 'SecHalfSession', 'Sparsity', 'InformationContent');
    [filepath, filename] = fileparts(ReadData{i});
    save([filepath filesep filename '_pathProperties.mat']);
    print(figure (1), '-djpeg', filename);
    close(figure (1));
    print(figure (2), '-djpeg', filename);
    close(figure (2));
    print(figure (3), '-djpeg', filename);
    close(figure (3));
    print(figure (4), '-djpeg', filename);
    close(figure (4));
end
