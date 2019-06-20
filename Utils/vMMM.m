%vMMM LB May 2019

%Purpose of this script is to evaluate one-dimensional circular data to
%estimate distribution parameters for multimodal or unimodal directional
%data. 

%Dependencies: 
%       - mvmdist: https://github.com/chrschy/mvmdist 
%       - CircStat2012a


% ADD PATHS 

%LOAD DATA
data=load('F:\ClarkP30_Recordings\ProcessedData\LB04_S20170603182951',...
    'events','frames','spikesID','Spikes','samplerate','ratemap','maze_size_cm');

for ii=1:size(data.events,2)
    
    for i=1:length(data.Spikes)
        
        disp([extractBefore(data.spikesID.TetrodeNum{i},'.'),...
            ' Cell ',num2str(data.spikesID.CellNum(i))])
        
        [data_video_spk,~]=createframes_w_spikebinary(temp,ii,i);
        spks_VEL=data_video_spk(data_video_spk(:,6)==1,4);
        
        [~,~,~,~,~,hdTuning(i,:,ii)]=tuningcurve(data_video_spk(data_video_spk(:,6)==0,4),spks_VEL,temp.samplerate);
        
        
    end
end

%RUN K-MEANS TO INITIALIZE PARAMETERS 
[idx, centers] = ckmeans(angles, nClusters);

mu1 = [1 2];
Sigma1 = [2 0; 0 0.5];
mu2 = [-3 -5];
Sigma2 = [1 0;0 1];
rng(1); % For reproducibility
X = [mvnrnd(mu1,Sigma1,1000);mvnrnd(mu2,Sigma2,1000)];

%ESTIMATE PARAMETERS BASED ON THE DATA 
p = [0.5; 0.5];         % Mixture weights.
mu = [-pi/2; pi/2];     % Component means.
kappa = [5; 10];        % Concentration parameters of components.

vmm = VonMisesMixture(p, mu, kappa);


%RUN EM 
