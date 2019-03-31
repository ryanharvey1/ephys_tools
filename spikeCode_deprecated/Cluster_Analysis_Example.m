%% Cluster_Analysis_Example
% Clustering Fisher's Iris Data Using K-Means Clustering

rng(14,'twister');
load fisheriris
%%
[cidx2,cmeans2] = kmeans(meas,2,'dist','sqeuclidean');
[silh2,h] = silhouette(meas,cidx2,'sqeuclidean');
%%
ptsymb = {'bs','r^','md','go','c+'};
for i = 1:2
    clust = find(cidx2==i);
    plot3(meas(clust,1),meas(clust,2),meas(clust,3),ptsymb{i});
    hold on
end
plot3(cmeans2(:,1),cmeans2(:,2),cmeans2(:,3),'ko');
plot3(cmeans2(:,1),cmeans2(:,2),cmeans2(:,3),'kx');
hold off
xlabel('Sepal Length');
ylabel('Sepal Width');
zlabel('Petal Length');
view(-137,10);
grid on
%%
[cidx3,cmeans3] = kmeans(meas,3,'display','iter');
%%
[cidx3,cmeans3,sumd3] = kmeans(meas,3,'replicates',5,'display','final');
%%
sum(sumd3)
%%
[silh3,h] = silhouette(meas,cidx3,'sqeuclidean');
%%
for i = 1:3
    clust = find(cidx3==i);
    plot3(meas(clust,1),meas(clust,2),meas(clust,3),ptsymb{i});
    hold on
end
plot3(cmeans3(:,1),cmeans3(:,2),cmeans3(:,3),'ko');
plot3(cmeans3(:,1),cmeans3(:,2),cmeans3(:,3),'kx');
hold off
xlabel('Sepal Length');
ylabel('Sepal Width');
zlabel('Petal Length');
view(-137,10);
grid on
%%
[mean(silh2) mean(silh3)]
%%
[cidxCos,cmeansCos] = kmeans(meas,3,'dist','cos');
%%
[silhCos,h] = silhouette(meas,cidxCos,'cos');
[mean(silh2) mean(silh3) mean(silhCos)]
%%
for i = 1:3
    clust = find(cidxCos==i);
    plot3(meas(clust,1),meas(clust,2),meas(clust,3),ptsymb{i});
    hold on
end
hold off
xlabel('Sepal Length');
ylabel('Sepal Width');
zlabel('Petal Length');
view(-137,10);
grid on
%%
lnsymb = {'b-','r-','m-'};
names = {'SL','SW','PL','PW'};
meas0 = meas ./ repmat(sqrt(sum(meas.^2,2)),1,4);
ymin = min(min(meas0));
ymax = max(max(meas0));
for i = 1:3
    subplot(1,3,i);
    plot(meas0(cidxCos==i,:)',lnsymb{i});
    hold on;
    plot(cmeansCos(i,:)','k-','LineWidth',2);
    hold off;
    title(sprintf('Cluster %d',i));
    xlim([.9, 4.1]);
    ylim([ymin, ymax]);
    h_gca = gca;
    h_gca.XTick = 1:4;
    h_gca.XTickLabel = names;
end

%%
subplot(1,1,1);
for i = 1:3
    clust = find(cidxCos==i);
    plot3(meas(clust,1),meas(clust,2),meas(clust,3),ptsymb{i});
    hold on
end
xlabel('Sepal Length');
ylabel('Sepal Width');
zlabel('Petal Length');
view(-137,10);
grid on
sidx = grp2idx(species);
miss = find(cidxCos ~= sidx);
plot3(meas(miss,1),meas(miss,2),meas(miss,3),'k*');
legend({'setosa','versicolor','virginica'});
hold off
%% Clustering Fisher's Iris Data Using Hierarchical Clustering

eucD = pdist(meas,'euclidean');
clustTreeEuc = linkage(eucD,'average');
%%
cophenet(clustTreeEuc,eucD)
%%
[h,nodes] = dendrogram(clustTreeEuc,0);
h_gca = gca;
h_gca.TickDir = 'out';
h_gca.TickLength = [.002 0];
h_gca.XTickLabel = [];
%%
cosD = pdist(meas,'cosine');
clustTreeCos = linkage(cosD,'average');
cophenet(clustTreeCos,cosD)
%%
[h,nodes] = dendrogram(clustTreeCos,0);
h_gca = gca;
h_gca.TickDir = 'out';
h_gca.TickLength = [.002 0];
h_gca.XTickLabel = [];
%%
[h,nodes] = dendrogram(clustTreeCos,12);
%%
[sum(ismember(nodes,[11 12 9 10])) sum(ismember(nodes,[6 7 8])) ...
    sum(ismember(nodes,[1 2 4 3])) sum(nodes==5)]
%%
hidx = cluster(clustTreeCos,'criterion','distance','cutoff',.006);
for i = 1:5
    clust = find(hidx==i);
    plot3(meas(clust,1),meas(clust,2),meas(clust,3),ptsymb{i});
    hold on
end
hold off
xlabel('Sepal Length');
ylabel('Sepal Width');
zlabel('Petal Length');
view(-137,10);
grid on

%%
clustTreeSng = linkage(eucD,'single');
[h,nodes] = dendrogram(clustTreeSng,0);
h_gca = gca;
h_gca.TickDir = 'out';
h_gca.TickLength = [.002 0];
h_gca.XTickLabel = [];
%%



