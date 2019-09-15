%Test cells 
%TT4, cell 3 (zero-order), 7 (trimodal), 4(unimodal); TT8 cell 3(unimodal)

%Load Data
data=load('F:\ClarkP30_Recordings\ProcessedData\LB04_S20170420130220');

%%
cellid={'TT4.mat',7};

cells_to_find=strcat(cellid{1},num2str(cellid{2}));

cell_list=strcat(data.spikesID.TetrodeNum,num2str(data.spikesID.CellNum));

cells=find(ismember(strrep(cell_list,' ',''),cells_to_find))';

%Pull Tuning Curves 
[data_video_spk,~]=createframes_w_spikebinary(data,1,cells);

[r,~,Ispk,~,pd,tuning]=tuningcurve(data_video_spk(data_video_spk(:,6)==0,4),...
    data_video_spk(data_video_spk(:,6)==1,4),data.samplerate);
fig = figure; 
fig.Color = [1 1 1];
postprocessFigures.plot_HD_tuning(data,1,cells)
title(sprintf('r: %4.2f DIC: %4.2f PD: %4.2f' ,[r,Ispk,pd]))
% export_fig(['D:\Users\BClarkLab\Google Drive (lberkowitz@unm.edu)\UNM PhD\Presentations\TalkFigs\',num2str(cells),'_CurrentApproach.png'],'-m4') 
%Initialize Parameters 

fig=figure; 
fig.Color=[1 1 1];
postprocessFigures.raster(data)
title(sprintf('r: %4.2f DIC: %4.2f PD: %4.2f' ,[r,Ispk,pd]))
% export_fig(['D:\Users\BClarkLab\Google Drive (lberkowitz@unm.edu)\UNM PhD\Presentations\TalkFigs\',num2str(cells),'_Raster.png'],'-m4') 

%Initialize Parameters 

% In this step, we need to create the mixture weights (p), spread of the
% distribution (kappa), and mean angles(mu). These variables are
% used as inputs into the VonMiseMixture.m function from the mvmdist-master
% toolbox. 

% Specify the max number of distributions you would like to test. 
 numDist = 8;

% Create mixture weights (p) whereby the summed weight of each distribution
% equals 1.

p{1} = 1; %initialize weight equal to 1
for i = 1:numDist
p{i+1} = repmat(1/i,i,1); %populate mixture weights as equal amounts adding to 1. 
end 

% Create the circular mean of each distribution (mu{1 x numDist}) in radians
% for each distribution set. 

mu = component_mu(numDist);

% Create kappa inputs
kappa{1} = 0; %initialize kappa for zero order model
for i = 1:numDist
kappa{i+1} = repmat(i,i,1); %populate mixture weights as equal amounts adding to 1. 
end

% 

vmm = VonMisesMixture(p{4}, mu{3}, kappa{4}); %creates our model 

samples=wrapToPi(deg2rad(data_video_spk(data_video_spk(:,6)==1,4)));
samples(isnan(samples))=[];

angles = linspace(-pi, pi, 1000)'; % The vmm.pdf() function expects a column-vector as input.
likelihoods = vmm.pdf(samples);

fig=figure;
fig.Color=[1 1 1];
scatter(samples, likelihoods,'k','filled'); grid on;
title('Data plotted as a function of the von Mises PDF')
% export_fig(['D:\Users\BClarkLab\Google Drive (lberkowitz@unm.edu)\UNM PhD\Presentations\TalkFigs\',num2str(cells),'_vMPdf.png'],'-m4') 

% Fit new model on data samples assuming n mixture components.
for i=1:8
fittedVmm = fitmvmdist(samples, i, ...
  'MaxIter', 250,'Replicates',5); % Set maximum number of EM iterations to 250 and number of parameter fitting replications to 5
fit_measures(i)=fittedVmm.logLikelihood;
end

% Plot initial and fitted distributions.
angles = linspace(-pi, pi, 1000)';
likelihoodsInitial = vmm.pdf(angles);
likelihoodsFitted = fittedVmm.pdf(angles);

%Likelihood figure 
fig=figure; 
fig.Color=[1 1 1];
% plot(angles, likelihoodsInitial,'k','LineWidth',2); hold on;
plot(angles, likelihoodsFitted,'r','LineWidth',2); grid on;
title('Likelihood values for initial and fitted models')
xlabel('Radians')
ylabel('LogLikeliihood Values')
axis([-pi, pi, 0, 1]);
export_fig(['D:\Users\BClarkLab\Google Drive (lberkowitz@unm.edu)\UNM PhD\Presentations\TalkFigs\',num2str(cells),'_Models.png'],'-m4') 


%Fitted Likelihoods over tuning curve
fig=figure;
fig.Color=[1 1 1]
postprocessFigures.plot_HD_tuning(data,1,cells)
hold on; 
Polarplot=polar(angles, rescale(likelihoodsFitted,0,max(tuning)));
set(Polarplot,'linewidth',2,'color','r');
title(['Preferred Directions: ',num2str(wrapTo360(rad2deg(fittedVmm.mu(1)))),', '...
    num2str(wrapTo360(rad2deg(fittedVmm.mu(2)))), ', ',num2str(wrapTo360(rad2deg(fittedVmm.mu(3))))])
export_fig(['D:\Users\BClarkLab\Google Drive (lberkowitz@unm.edu)\UNM PhD\Presentations\TalkFigs\',num2str(cells),'_NewApproach.png'],'-m4') 



