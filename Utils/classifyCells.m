% classifyCells
%
%Exploratory classification analysis to discern different functional cell
%types recorded from posterior cortex. Created for Multivariate Statistics
%psych 650,  LB 2019

data=compileResults('F:\ClarkP30_Recordings\ProcessedData');

group={'LB03','LB05','LB04','LB06','LB07'};

%% COMPILE GROUP ID FOR ALL SESSIONS
data.measures=[];
data.id=[];
for i=1:length(group)
    data.measures=cat(1,data.measures,data.(group{i}).measures);
    data.id=cat(1,data.id,data.(group{i}).id);
end

% COMPILE DATA
group=data.measures(:,:,1);
groupid=data.id;


%% DELETE MEASURES FOR LINEAR TRACK
varnames=data.varnames(:,[1:11,30,37,43:46,60]); group=group(:,[1:11,30,37,43:46,60]);

%COMPUTE COMMON HD MEASURES
HDvar={'Stability'};
for i=1:length(groupid)
    
    temp=load(['F:\ClarkP30_Recordings\ProcessedData\',groupid{i,1}],...
        'events','frames','spikesID','Spikes','samplerate','ratemap','maze_size_cm','thetaautocorr');
    
    cell=find(contains(temp.spikesID.TetrodeNum,groupid{i,2}) & ismember(temp.spikesID.CellNum,str2double(groupid{i,3})))';
    
    autocors(i,:)=temp.thetaautocorr{cell,1};
    
    [data_video_spk,~]=createframes_w_spikebinary(temp,1,cell); %get spike data from first session only
    spks_VEL=data_video_spk(data_video_spk(:,6)==1,4);
    
    %meas varnames = rlength,information/spk,stability,distributive
    %ratio.
    
    [~,~,~,~,~,hdTuning(i,:)]=tuningcurve(data_video_spk(data_video_spk(:,6)==0,4),spks_VEL,temp.samplerate);
    
    meas(i,1)=HD_cell_analysis.hd_stability(data_video_spk);
               
end

%% Center Tuning curves to evaluate shape 

for i=1:size(hdTuning,1)

    [M,I]=max(hdTuning(i,:));
    middlebin=round(median(1:length(hdTuning(i,:))));
    
    hdTuningshift=circshift(hdTuning(i,:),(middlebin-I)-1);
    normTuning(i,:)=rescale(hdTuningshift,0,1);
    
end

clear middlebin M I hdTuningshift

%% Perform PCA to examine variables 
var=[varnames HDvar];
meas=[group meas];

% normalize data matrix column-wise to get variables on same scale. 
for i=1:size(meas,2)
normMeas(:,i)=rescale(meas(:,i));
end
sum(isnan(normMeas),2)>0

test=normMeas;

test(sum(isnan(normMeas),2)>0,:)=[]; %remove missing data

[coeff,score,latent,tsquared,explained,mu] = pca(test); %run PCA on normalized data 

fig=figure; 
biplot(coeff(:,1:2),'Scores',score(:,1:2),'Marker','o','MarkerFaceColor','r','VarLabels',var);% ,'VarLabels',var
darkBackground(fig,[0.2 0.2 0.2],[1 1 1])

fig=figure; fig.Color=[1 1 1]; 
plot(explained,'LineWidth',2,'Color','r')
title('Scree Plot')
ylabel('Explained Variance (%)')
xlabel('Components')

darkBackground(fig,[0.2 0.2 0.2],[1 1 1])

fig=figure; 
scatter3(score(:,1),score(:,2),score(:,3),'MarkerFaceColor','r','MarkerEdgeColor','r','MarkerFaceAlpha',.5)
axis equal
xlabel('1st Principal Component')
ylabel('2nd Principal Component')
zlabel('3rd Principal Component')%% Clustering Analysis 
darkBackground(fig,[0.2 0.2 0.2],[1 1 1])


%% Clustering the data

clust = zeros(size(test,1),10);
for i=1:10
clust(:,i) = kmeans(test,i,'emptyaction','singleton',...
        'replicate',5);
end

%EVALUATE K_MEANS WITH SILHOUETTE 
eva = evalclusters(test,clust,'silhouette') 

fig=figure; 
silhouette(test,clust(:,2))

fig=figure; %LABELED COMPONENTS
scatter3(score(clust(:,2)==1,1),score(clust(:,2)==1,2),score(clust(:,2)==1,3),'MarkerFaceColor','r','MarkerEdgeColor','r','MarkerFaceAlpha',.5)
hold on
scatter3(score(clust(:,2)==2,1),score(clust(:,2)==2,2),score(clust(:,2)==2,3),'MarkerFaceColor','b','MarkerEdgeColor','b','MarkerFaceAlpha',.5)
hold off
axis equal
xlabel('1st Principal Component')
ylabel('2nd Principal Component')
zlabel('3rd Principal Component')%% Clustering Analysis 
darkBackground(fig,[0.2 0.2 0.2],[1 1 1])

fig=figure; %AVERAGE NORMALIZED TUNING CURVE BETWEEN CLUSTERS
subplot(1,2,1); 
plot(nanmean(testTuning(clust(:,2)==1,:))','LineWidth',2,'Color','r'); title('Cluster 1: Average Tuning Curve');
subplot(1,2,2)
plot(nanmean(testTuning(clust(:,2)==2,:))','LineWidth',2,'Color','b');title('Cluster 1: Average Tuning Curve');
darkBackground(fig,[0.2 0.2 0.2],[1 1 1])

corr2(nanmean(testTuning(clust(:,2)==1,:)),nanmean(testTuning(clust(:,2)==2,:)));
fig=figure; %POP VECTORS BETWEEN CLUSTERS
subplot(1,2,1); 
imagesc(testTuning(clust(:,2)==1,:)); title('Cluster 1: Population Vectors');
subplot(1,2,2)
imagesc(testTuning(clust(:,2)==2,:));title('Cluster 1: Population Vectors');
darkBackground(fig,[0.2 0.2 0.2],[1 1 1])
% fig=figure; %AVERAGE NORMALIZED TUNING CURVE BETWEEN CLUSTERS
% subplot(1,2,1); 
% plot(nanmean(normTuning(clust(:,2)==1,:))','LineWidth',2,'Color','r'); title('Cluster 1: Average Tuning Curve');
% subplot(1,2,2)
% plot(nanmean(normTuning(clust(:,2)==2,:))','LineWidth',2,'Color','b');title('Cluster 1: Average Tuning Curve');
% darkBackground(fig,[0.2 0.2 0.2],[1 1 1])
% 
% fig=figure; %POP VECTORS BETWEEN CLUSTERS
% subplot(1,2,1); 
% imagesc(normTuning(clust(:,2)==1,:)); title('Cluster 1: Population Vectors');
% subplot(1,2,2)
% imagesc(normTuning(clust(:,2)==2,:));title('Cluster 1: Population Vectors');
% darkBackground(fig,[0.2 0.2 0.2],[1 1 1])


figure; %COMPARE MEASURES
testMeas=meas;
testMeas(sum(isnan(normMeas),2)>0,:)=[];
CDFplots(testMeas(clust(:,2)==1,:),testMeas(clust(:,2)==2,:),{'cluster 1','cluster 2'},vars,1)

%% CLUSTER DATA WITH GUASSIAN MIXTURE MODELLING With REGULARIZATION -- choose GMM with lowest AIC

%%
rng(3)
k = 2;
d = 500;
x1 = linspace(min(score(:,1)) - 2,max(score(:,1)) + 2,d);
x2 = linspace(min(score(:,2)) - 2,max(score(:,2)) + 2,d);
[x1grid,x2grid] = meshgrid(x1,x2);
X0 = [x1grid(:) x2grid(:)];
threshold = sqrt(chi2inv(0.99,2));
options = statset('MaxIter',1000); % Increase number of EM iterations
AIC = zeros(2,1);
gmfit = cell(2,1);
figure;

for k=1:2
gmfit{k} = fitgmdist(score(:,1:2),k,'CovarianceType','diagonal',...
    'SharedCovariance',true,'Options',options);
clusterX = cluster(gmfit{k},score(:,1:2));
mahalDist = mahal(gmfit{k},X0);
subplot(1,2,k);
h1 = gscatter(score(:,1),score(:,2),clusterX);
hold on;
for m = 1:k
    idx = mahalDist(:,m)<=threshold;
    Color = h1(m).Color*0.75 + -0.5*(h1(m).Color - 1);
    h2 = plot(X0(idx,1),X0(idx,2),'.','Color',Color,'MarkerSize',1);
    uistack(h2,'bottom');
end
plot(gmfit{k}.mu(:,1),gmfit{k}.mu(:,2),'kx','LineWidth',2,'MarkerSize',10)
title(sprintf('Sigma is %s, SharedCovariance = %s',...
    'diagonal',true),'FontSize',8)
legend(h1,{'1','2','3'});
hold off
AIC(k)= GMModels{k}.AIC;
end
darkBackground(gcf,[0.2 0.2 0.2],[1 1 1])

[minAIC,numComponents] = min(AIC);
numComponents

FisrtModel = GMModels{1}