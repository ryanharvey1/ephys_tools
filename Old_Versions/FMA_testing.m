addpath('/Users/ryanharvey/GoogleDrive/MatlabDir/BClarkToolbox/Analysis');
path='/Users/RyanHarvey/GoogleDrive/postprocessMclust_test/RH13/2016-07-02_18-35-24';
% path='/Volumes/SamsungUSB/TESTSESSION/VELFILTER/2017-05-16_15-16-32';
% path='/Users/ryanharvey/GoogleDrive/Takeout/2017-04-16_18-05-19';
% path='F:\Users\reharvey\Place_Cell_Data\PAE_Rat\LS23\2017-04-16_18-05-19';
postprocessMClust_v9working2(path,120,'yes',true)
%%

% ratemap creation Linear track

x=(data_video_smoothfilt(:,2)-min(data_video_smoothfilt(:,2)))/range(data_video_smoothfilt(:,2));
y=(data_video_smoothfilt(:,3)-min(data_video_smoothfilt(:,3)))/range(data_video_smoothfilt(:,3));

nBinsx = round(track_length/5); nBinsy = 1;

map = Map([data_video_smoothfilt(:,1) x],spks_VEL(:,1),'nBins',[nBinsx nBinsy]);

stats = MapStats(map);

% [map,stats] = FiringMap([data_video_smoothfilt(:,1) x y],spks_VEL(:,1));

 figure;
 figure(1);PlotColorMap(map.z,map.time); 
 figure(2);pcolor([map.z;map.z]); shading interp; colormap jet;
 %% Arena 
x=(data_video_smoothfilt(:,2)-min(data_video_smoothfilt(:,2)))/range(data_video_smoothfilt(:,2));
y=(data_video_smoothfilt(:,3)-min(data_video_smoothfilt(:,3)))/range(data_video_smoothfilt(:,3));

nBinsx = round(track_length/5); nBinsy = nBinsx;

map = Map([data_video_smoothfilt(:,1) x y],spks_VEL(:,1),'nBins',[nBinsx nBinsy]);
matrix=zeros(nBinsx+1, nBinsy+1); workingmap=map.z; matrix(1:nBinsx,1:nBinsy)=workingmap;
[ circBinMap ] = circRateMap(matrix,nBinsx+1 );
circBinMap(end,:) = [];
circBinMap(:,end) = [];
map.z=circBinMap;
stats = MapStats(map);

% [map,stats] = FiringMap([data_video_smoothfilt(:,1) x y],spks_VEL(:,1));

 figure;
 figure(1);PlotColorMap(map.z,map.time); 
 figure(2);pcolor(map.z); shading interp; colormap jet;
 figure(3);pcolor(circBinMap); shading interp; colormap jet;
 
 %%
 % Test velocity filtering 
 x=(data_video_smoothfilt(:,2)-min(data_video_smoothfilt(:,2)))/range(data_video_smoothfilt(:,2));
y=(data_video_smoothfilt(:,3)-min(data_video_smoothfilt(:,3)))/range(data_video_smoothfilt(:,3));

 V = LinearVelocity([[1:length(data_video_smoothfilt(:,1))]' x y]);
 
 [periods,quiescence] = QuietPeriods([data_video_smoothfilt(:,1) V*30],5,3);