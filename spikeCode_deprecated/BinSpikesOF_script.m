% function BinSpikesOF
% loads tetrodes clusters, restricts them to run times

%%
% load variables from file (even if computed from scratch and saved above)
load runpos;
%%
dispFig=input('do you want to see plotted figures?','s');
%% change to data directory
if ~exist('strawberryfields', 'dir')
	fprintf('cannot find hard-coded data directory - strawberryfields\n');
	fprintf('creating directory for you - populate it with .t files please\n');
	mkdir('strawberryfields');
	cd(startdir);
	return;
end
cd('strawberryfields');

% Bin the positional data into ~50 (a couple extra so we can trim the sides 
% as necessary bins along the X and Y axes
nBins=25; % we will cut off 3 bins on either side to deal with areas not really travelled to

MinY=min(min(Data(Y_run1)),min(Data(Y_run2)));
MaxY=max(max(Data(Y_run1)),max(Data(Y_run2)));

MinX=min(min(Data(X_run1)),min(Data(X_run2)));
MaxX=max(max(Data(X_run1)),max(Data(X_run2)));

epsilon = 0.0001;

edges{1}=linspace(MinY,MaxY+epsilon,nBins+1);
edges{2}=linspace(MinX,MaxX+epsilon,nBins+1);

dat4matrix1=[Data(Y_run1), Data(X_run1)];
dat4matrix2=[Data(Y_run2), Data(X_run2)];

% bin positional data
matrix1=flipud(hist3(dat4matrix1,'Edges',edges));
matrix1(1,:)=[];
matrix1(:,end)=[];
matrix2=flipud(hist3(dat4matrix2,'Edges',edges));
matrix2(1,:)=[];
matrix2(:,end)=[];

% Convert counts to occupancy
dt1=DT(X1)/1e4; % calculate frame interval and convert to sec from 100 us units 
dt2=DT(X2)/1e4;

Occ1=matrix1*dt1; % occupancy (in sec) equals counts times frame interval 
Occ2=matrix2*dt2;

% find pixels with zero occupancy- these will get blanked out out in the
% rate maps
Occ1NonZeroLogical=Occ1~=0;
Occ2NonZeroLogical=Occ2~=0;
Occ1ZeroLogical=Occ1==0;
Occ2ZeroLogical=Occ2==0;
 
% plot path of rat
h=figure;

subplot(2,2,1);
scatter(dat4matrix1(:,2),dat4matrix1(:,1),'.')
title('path run 1')

subplot(2,2,2);
scatter(dat4matrix2(:,2),dat4matrix2(:,1),'.')
title('path run 2')

% make occupancy heat map
subplot(2,2,3);
imagesc(Occ1)
title('occupancy run 1')

subplot(2,2,4);
imagesc(Occ2)
title('occupancy run 2')

saveas(h,'path and occupancy');
close

% Load spikes
tfiles = FindFiles('*.t');  % looks in folder you are in and all subfolders too
S = LoadSpikes(tfiles);
nCells=length(S);

cellID=zeros(nCells, 2);  % array with the ID of each cell. Tetrode number col 1, cell number col 2.
S_maze1=cell(nCells, 1);
S_maze2=cell(nCells, 1);
S_run1=cell(nCells, 1);
S_run2=cell(nCells, 1);

% create arrays for rate maps (i.e. raw, filled, and smoothed, all for run1 and run2)
RawRateMaps1=zeros(nBins,nBins,nCells);
RawRateMaps2=zeros(nBins,nBins,nCells);

FilledRateMaps1=zeros(nBins,nBins,nCells);
FilledRateMaps2=zeros(nBins,nBins,nCells);

SmoothedRateMaps1=zeros(nBins,nBins,nCells);
SmoothedRateMaps2=zeros(nBins,nBins,nCells);

for c=1:nCells
	% restrict S to just maze epoch
	S_maze1{c} = Restrict(S{c}, startmaze, endmaze);  % Restrict cell spikes to lap epoch time range
	S_run1{c} = Restrict(S_maze1{c}, startrun1, stoprun1);
	clear S_maze1;
	S_maze2{c} = Restrict(S{c}, startmaze, endmaze);  
	S_run2{c} = Restrict(S_maze2{c}, startrun2, stoprun2);
	clear S_maze2;
	% find tetrode and cell number and store in cellID array
	[pathstr,name,ext]=fileparts(tfiles{c});
	cellID(c,1)=str2num(name(3:4));
    if name(5)=='d'
        cellID(c,2)=str2num(name(7:end));
    elseif name(5)=='c'
        cellID(c,2)=str2num(name(9:end));
    end
    
    % interpolate the position of spikes
    X1_spikes=interp1(Range(X_run1), Data(X_run1), Range(S_run1{c}), 'linear'); % Find location of spikes
    Y1_spikes=interp1(Range(Y_run1), Data(Y_run1), Range(S_run1{c}), 'linear');
    
    X2_spikes=interp1(Range(X_run2), Data(X_run2), Range(S_run2{c}), 'linear');
    Y2_spikes=interp1(Range(Y_run2), Data(Y_run2), Range(S_run2{c}), 'linear'); 
    
    dat4Smatrix1=[Y1_spikes, X1_spikes];
    dat4Smatrix2=[Y2_spikes, X2_spikes];
    
    % make matrix of binned spikes
    
    Smatrix1=flipud(hist3(dat4Smatrix1,'Edges',edges));
    Smatrix1(1,:)=[];
    Smatrix1(:,end)=[];
    Smatrix2=flipud(hist3(dat4Smatrix2,'Edges',edges));
    Smatrix2(1,:)=[];
    Smatrix2(:,end)=[];
    
    % divide binned spikes by occupancy to get rate maps and store data
    BFiringRateMatrix1=Smatrix1./Occ1;
    BFiringRateMatrix2=Smatrix2./Occ2; 
    
    RawRateMaps1(:,:,c)=BFiringRateMatrix1;
    RawRateMaps2(:,:,c)=BFiringRateMatrix2;
    
    % Fill in NaN's (bins where Occ was ==0 and no spikes fired), except outside of track
    % edges
    FilledFRMatrix1=FillNaNs(BFiringRateMatrix1);
    FilledFRMatrix2=FillNaNs(BFiringRateMatrix2);
    
    % Fill in Inf's (bins where Occ was ==0, but spikes were assigned due to
    % discontinuities in video tracking and interpolation of spike
    % locations
    
    FilledFRMatrix1=FillInfs(FilledFRMatrix1);
    FilledFRMatrix2=FillInfs(FilledFRMatrix2);
    
    % Input values in data array
    FilledRateMaps1(:,:,c)=FilledFRMatrix1;
    FilledRateMaps2(:,:,c)=FilledFRMatrix2;
    
    % Smooth rate map with 5x5 boxcar filter

    SmoothedRateMap1=Smooth2D(FilledFRMatrix1);
    SmoothedRateMap2=Smooth2D(FilledFRMatrix2);
    
    SmoothedRateMaps1(:,:,c)=SmoothedRateMap1;
    SmoothedRateMaps2(:,:,c)=SmoothedRateMap2;
    
    % Normalize the smoothed rate map to the median of the highest nine
    % values (from Leutgeb paper)- same as taking the 5th highest value
    SRMNzero1=SmoothedRateMap1(Occ1NonZeroLogical);
    SRMNzero1Ind=find(Occ1NonZeroLogical);
    [RateMap1Vals,RateMap1Inds]=sort(SRMNzero1,'descend');
    RateMap1Max=RateMap1Vals(5);
    SmoothedRateMap1(SRMNzero1Ind(RateMap1Inds(1:4)))=RateMap1Max;    

    SRMNzero2=SmoothedRateMap2(Occ2NonZeroLogical);
    SRMNzero2Ind=find(Occ2NonZeroLogical);
    [RateMap2Vals,RateMap2Inds]=sort(SRMNzero2,'descend');
    RateMap2Max=RateMap2Vals(5);
    SmoothedRateMap2(SRMNzero2Ind(RateMap2Inds(1:4)))=RateMap2Max;    

    % scale data so that occupied pixels are imaged with jet cmap and
    % unoccupied pixels are white. This is accomplished by making values 
    % from occupied pixels betweeen 1 and 64, and unoccupied pixels =65.
    % code is adapted from http://www.mathworks.com/support/tech-notes/1200/1215.html
    
    m = 64;  % 64-elements is each colormap
    
    % scale data from occupied bins to the first colour map, from 
    % unoccupied bins to the second map. Map will be: colormap([jet(64);bone(64)])

    ScaledSRM1=zeros(nBins,nBins);
    ScaledSRM2=zeros(nBins,nBins);   
    
    % the 2 rate maps will be scaled the same (i.e. using same max and min
    % values so that they can be directly compared)
    RateMap1Min=min(SRMNzero1);
    RateMap2Min=min(SRMNzero2);
    RateMapMin=min(RateMap1Min,RateMap2Min);
    RateMapMax=max(RateMap1Max,RateMap2Max);
    
    ScaledSRM1(Occ1NonZeroLogical) = min(m,round((m-1)*(SmoothedRateMap1(Occ1NonZeroLogical)-RateMapMin)/(RateMapMax-RateMapMin))+1);
    ScaledSRM2(Occ2NonZeroLogical) = min(m,round((m-1)*(SmoothedRateMap2(Occ2NonZeroLogical)-RateMapMin)/(RateMapMax-RateMapMin))+1);
    
    ScaledSRM1(Occ1ZeroLogical)=128; % will be white (i.e. last value in bone colormap)
    ScaledSRM2(Occ2ZeroLogical)=128;  
    
    % create figure
    htitle=['TT ', num2str(cellID(c,1)), ' c ', num2str(cellID(c,2))];
    if dispFig=='y'
            h=figure('Name',htitle,'NumberTitle','off', 'Units', 'pixels', 'Position', [0 0 1000 500]);    
    else
        h=figure('Name',htitle,'NumberTitle','off', 'Units', 'pixels', 'Position', [0 0 1000 500], 'Visible', 'off');
    end
    set(h, 'Unit', 'pixels')
    set(h,'position', [50, 50, 1200, 600])
    % set(h,'Renderer','OpenGL')
    hold on
    
    % plot spikes over path
    j=subplot(1,2,1);
    set(j, 'Unit', 'pixels')
    set(j,'position', [50, 50, 400, 400])
    plot(dat4matrix1(:,2),dat4matrix1(:,1),'*b', 'MarkerSize',2)
    axis([MinX, MaxX, MinY, MaxY])   
    hold on
    % scatter(X1_spikes,Y1_spikes,10,'r', 'filled')
    plot(X1_spikes,Y1_spikes,'.r', 'MarkerSize',8)
    title('spikes run 1')

    j=subplot(1,2,2);
    set(j, 'Unit', 'pixels')
    set(j,'position', [550, 50, 400, 400])
    plot(dat4matrix2(:,2),dat4matrix2(:,1),'*b', 'MarkerSize',2)
    axis([MinX, MaxX, MinY, MaxY])   
    hold on
    % scatter(X2_spikes,Y2_spikes,10,'r', 'filled')
    plot(X2_spikes,Y2_spikes,'.r', 'MarkerSize',8)
    title('spikes run 2')
    
    fstrg = ['OpenField_TT', num2str(cellID(c,1)), 'c', num2str(cellID(c,2)),'_scatter.tiff']; 
	% print(h,'-dtiff', fstrg)    
    saveSameSize(h,'format','tiff', 'file', fstrg);
    close   
    
    % create rate map figure
    htitle=['TT ', num2str(cellID(c,1)), ' c ', num2str(cellID(c,2))];
    if dispFig=='y'
        h=figure('Name',htitle,'NumberTitle','off','Units', 'pixels', 'Position', [0 0 1050 500]);
    else
                h=figure('Name',htitle,'NumberTitle','off','Units', 'pixels', 'Position', [0 0 1050 500],'Visible', 'off');
    end
    hold on
    % make rate map tif
    j=subplot(1,2,1);
    set(j, 'Unit', 'pixels')
    set(j,'position', [50, 50, nBins*16,nBins*16])
    colormap([jet(64);bone(64)])
    image(1:16*(nBins), 1:16*(nBins), ScaledSRM1)
    title('rate map run 1')
    
    j=subplot(1,2,2);
    set(j, 'Unit', 'pixels')
    set(j,'position', [575, 50, nBins*16,nBins*16])
    colormap([jet(64);bone(64)])
    image(1:16*(nBins), 1:16*(nBins), ScaledSRM2)
    title('rate map run 1')
    
    fstrg = ['OpenField_' 'TT', num2str(cellID(c,1)), 'c', num2str(cellID(c,2)),'_RateMap']; 
	% print(h,'-dtiff', fstrg)    
    saveSameSize(h,'format','tiff', 'file', fstrg);
    close   
    
end
clear S 

save('RateMaps.mat', 'Occ1', 'Occ2', 'cellID', 'RawRateMaps1', 'RawRateMaps2',...
    'FilledRateMaps1', 'FilledRateMaps2', 'SmoothedRateMaps1',...
    'SmoothedRateMaps2', 'nBins');

clear




        
    
        
    





