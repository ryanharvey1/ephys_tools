% analyze_thetaautocorr
% this code loads in Ryan's post process data structure and compiles 
% theta autocorrelations

% load('data.mat')
close all
control={'RH13','RH14','LS21','LS23','LE2821','LE2823'};
PAE={'RH11','RH16','LS17','LS19','LE2813'};

% EXTRACT MEASURES FOR LINEAR TRACK
[group1,idc]=extractLINEAR(control,data);
[group2,idp]=extractLINEAR(PAE,data);

% EXTRACT MEASURES FOR CYLINDER
[group1cylinder]=extractCYLINDER(control,data);
[group2cylinder]=extractCYLINDER(PAE,data);

% FILTER OUT NON-PLACE CELLS
[group1,idc]=placefilter(group1,idc);
[group2,idp]=placefilter(group2,idp);
% [group1cylinder]=placefilter(group1cylinder);
% [group2cylinder]=placefilter(group2cylinder);

autoc=extractautocorr(idc,data);
autop=extractautocorr(idp,data);

% SORT AND PLOT MATRIX
[s,I]=sort(group1(:,35));
autoc=autoc(I,:);
autoc=autoc(~isnan(autoc(:,1)),:);


[s,I]=sort(group2(:,35));
autop=autop(I,:);
autop=autop(~isnan(autop(:,1)),:);


figure;imagesc(autoc);colormap jet;hold on;box off;axis off;
plot(rescale(-nansum(autoc),1,size(autoc,1)),'w','LineWidth',3)

figure;imagesc(autop);colormap jet;hold on;box off;axis off;
plot(rescale(-sum(autop),1,size(autop,1)),'w','LineWidth',3)


% \\\\\\\\\\\\\\\\\\\\ LOCAL FUNCTIONS BELOW \\\\\\\\\\\\\\\\\\\\
function [cells,ids]=extractLINEAR(rats,data)
cells=[];
ids=[];
for r=1:length(rats)
    sessions=fieldnames(data.(rats{r}));
    for s=1:length(sessions)
        ncells=fieldnames(data.(rats{r}).(sessions{s}).('thetaautocorr'));
        ids=[ids;repmat({(rats{r}),(sessions{s}),('thetaautocorr')},length(ncells),1),ncells,repmat({'session1'},length(ncells),1)];
        ids=[ids;repmat({(rats{r}),(sessions{s}),('thetaautocorr')},length(ncells),1),ncells,repmat({'session2'},length(ncells),1)];
        cells=[cells;[data.(rats{r}).(sessions{s}).measures(:,:,1);data.(rats{r}).(sessions{s}).measures(:,:,2)]];
    end
end
end

function [cells]=extractCYLINDER(rats,data)
cells=[];
for r=1:length(rats)
    sessions=fieldnames(data.(rats{r}));
    for s=1:length(sessions)
        [~,~,d]=size(data.(rats{r}).(sessions{s}).measures);
        if d>3 
            cells=[cells;data.(rats{r}).(sessions{s}).measures(:,:,[3:4])];
        end
    end
end
end

function [groupout,id]=placefilter(groupin,id)
groupout=groupin(groupin(:,1,1)>.3 & groupin(:,8,1)>50 & groupin(:,32,1)>groupin(:,8,1)*.2 & groupin(:,25,1)>0,:,:); 
id=id(groupin(:,1,1)>.3 & groupin(:,8,1)>50 & groupin(:,32,1)>groupin(:,8,1)*.2 & groupin(:,25,1)>0,:);

        cells=groupin(groupin(:,1,1:2)>0.25 &... % info content >.25
            data.(rats{r}).(sessions{s}).measures(:,5,direc)>0.3 &... % >.3 hz average rate
            data.(rats{r}).(sessions{s}).measures(:,4,direc)>2 &...
            data.(rats{r}).(sessions{s}).measures(:,8,direc)>50 &... % >50 spikes
            data.(rats{r}).(sessions{s}).measures(:,38,direc)<.05 &... % <0.5% spikes <2ms in autocorr
            data.(rats{r}).(sessions{s}).measures(:,39,direc)<0.5 &... % <50% cut off by threshold
            data.(rats{r}).(sessions{s}).measures(:,54,direc)>5 &...% above 5 laps
            data.(rats{r}).(sessions{s}).measures(:,53,direc)>.02); % temporal stability >.02

end

function auto=extractautocorr(id,data)
for i=1:length(id)
    if isfield(data.(id{i,1}).(id{i,2}).(id{i,3}).(id{i,4}),(id{i,5}))
        auto(i,:)=data.(id{i,1}).(id{i,2}).(id{i,3}).(id{i,4}).(id{i,5});
    end
end
end