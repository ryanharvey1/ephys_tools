% AnalyzeResults_HPCatn_lesion
data=compileResults('F:\Projects\HPCatn\ProcessedData');

control={'HPCatn02','HPCatn05'};
lesion={'HPCatn03','HPCatn04'};

% % load inactivation data (HPCatn06 / HPCatn07)
% load('F:\Projects\HPCatn\AnimalMetadata\HPCatn06_metadata.mat')
% clear mazes
% sessions=fieldnames(AnimalMetadata.RecordingLogs);
% for i=1:length(sessions)
%     mazes{i}=AnimalMetadata.RecordingLogs.(sessions{i}).MazeTypes;
% end
% idx=contains(data.HPCatn06.id(:,1),sessions(contains(mazes,'track')));
% 
% % baseline & 6 hour recovery
% HPCatn06_baseline=[data.HPCatn06.measures(idx,:,1);data.HPCatn06.measures(idx,:,2)];
% HPCatn06_baseline=[HPCatn06_baseline,[ones(size(data.HPCatn06.measures(idx,:,1),1),1);ones(size(data.HPCatn06.measures(idx,:,1),1),1)+1]];
% HPCatn06_baseline_id=[data.HPCatn06.id(idx,:);data.HPCatn06.id(idx,:)];
% 
% % inactivation
% HPCatn06_inactivation=[data.HPCatn06.measures(idx,:,3);data.HPCatn06.measures(idx,:,4)];
% HPCatn06_inactivation=[HPCatn06_inactivation,[ones(size(data.HPCatn06.measures(idx,:,1),1),1);ones(size(data.HPCatn06.measures(idx,:,1),1),1)+1]];
% HPCatn06_inactivation_id=[data.HPCatn06.id(idx,:);data.HPCatn06.id(idx,:)];
% 
% load('F:\Projects\HPCatn\AnimalMetadata\HPCatn07_metadata.mat')
% clear mazes
% sessions=fieldnames(AnimalMetadata.RecordingLogs);
% for i=1:length(sessions)
%     mazes{i}=AnimalMetadata.RecordingLogs.(sessions{i}).MazeTypes;
% end
% idx=contains(data.HPCatn07.id(:,1),sessions(contains(mazes,'track')));
% 
% % baseline & 6 hour recovery
% HPCatn07_baseline=[data.HPCatn07.measures(idx,:,1);data.HPCatn07.measures(idx,:,2)];
% HPCatn07_baseline=[HPCatn07_baseline,[ones(size(data.HPCatn07.measures(idx,:,1),1),1);ones(size(data.HPCatn07.measures(idx,:,1),1),1)+1]];
% HPCatn07_baseline_id=[data.HPCatn07.id(idx,:);data.HPCatn07.id(idx,:)];
% 
% % inactivation
% HPCatn07_inactivation=[data.HPCatn07.measures(idx,:,3);data.HPCatn07.measures(idx,:,4)];
% HPCatn07_inactivation=[HPCatn07_inactivation,[ones(size(data.HPCatn07.measures(idx,:,1),1),1);ones(size(data.HPCatn07.measures(idx,:,1),1),1)+1]];
% HPCatn07_inactivation_id=[data.HPCatn07.id(idx,:);data.HPCatn07.id(idx,:)];

%% COMPILE GROUPS
data.control.measures=[];
data.control.id=[];
for i=1:length(control)
    data.control.measures=cat(1,data.control.measures,data.(control{i}).measures);
    data.control.id=cat(1,data.control.id,data.(control{i}).id);
end

data.lesion.measures=[];
data.lesion.id=[];
for i=1:length(lesion)
    data.lesion.measures=cat(1,data.lesion.measures,data.(lesion{i}).measures);
    data.lesion.id=cat(1,data.lesion.id,data.(lesion{i}).id);
end
%% COMPILE LINEAR TRACK DATA
% also add index for running direction
group1=[[data.control.measures(:,:,1),ones(size(data.control.measures(:,:,1),1),1)];...
    [data.control.measures(:,:,2),ones(size(data.control.measures(:,:,1),1),1)+1]];

group2=[[data.lesion.measures(:,:,1),ones(size(data.lesion.measures(:,:,1),1),1)];...
    [data.lesion.measures(:,:,2),ones(size(data.lesion.measures(:,:,1),1),1)+1]];

group1id=[data.control.id;data.control.id];
group2id=[data.lesion.id;data.lesion.id];

% %% add baseline control data
% group1=[group1;HPCatn06_baseline;HPCatn07_baseline];
% group1id=[group1id;HPCatn06_baseline_id;HPCatn07_baseline_id];
% control{3}='HPCatn06';
% control{4}='HPCatn07';
% 
% %% add inactivation data to lesion group
% group2=[group2;HPCatn06_inactivation;HPCatn07_inactivation];
% group2id=[group2id;HPCatn06_inactivation_id;HPCatn07_inactivation_id];
% lesion{3}='HPCatn06';
% lesion{4}='HPCatn07';

%% DELETE MEASURES FOR OPEN ARENA
varnames=data.varnames;
varnames=[varnames,'runningdir'];
% colstodelete=size(group1,1)==sum(isnan(group1));
colstodelete=contains(varnames,["Cluster Grade","borderScore","E",...
    "DisplacementCorr","bordermod","egomod","Tightness","Incompleteness",...
    "StationInTime","TempMatch","BDistanceClust","BDistanceSpike"]);
varnames(colstodelete)=[];
group1(:,colstodelete,:)=[];
group2(:,colstodelete,:)=[];


%% SPLIT BY REGION
% load metadata files and extract region info
group1id=get_region_id(group1id,'F:\Projects\HPCatn\AnimalMetadata');
group1id(contains(group1id(:,end),'dg'),end) = cellstr(repmat('ca3',sum(contains(group1id(:,end),'dg')),1));

group2id=get_region_id(group2id,'F:\Projects\HPCatn\AnimalMetadata');
group2id(contains(group2id(:,end),'dg'),end) = cellstr(repmat('ca3',sum(contains(group2id(:,end),'dg')),1));



group1ca1 = group1(strcmp(group1id(:,4),'ca1'),:);
group1ca1id = group1id(strcmp(group1id(:,4),'ca1'),:);
group1ca3 = group1(strcmp(group1id(:,4),'ca3'),:);
group1ca3id = group1id(strcmp(group1id(:,4),'ca3'),:);


group2ca1 = group2(strcmp(group2id(:,4),'ca1'),:);
group2ca1id = group2id(strcmp(group2id(:,4),'ca1'),:);
group2ca3 = group2(strcmp(group2id(:,4),'ca3'),:);
group2ca3id = group2id(strcmp(group2id(:,4),'ca3'),:);



%% PLACE CELL FILTER
[group1ca1,group1ca1id]=placefieldfilter(group1ca1,group1ca1id,varnames);
[group2ca1,group2ca1id]=placefieldfilter(group2ca1,group2ca1id,varnames);
[group1ca3,group1ca3id]=placefieldfilter(group1ca3,group1ca3id,varnames);
[group2ca3,group2ca3id]=placefieldfilter(group2ca3,group2ca3id,varnames);



[uCA,~,~] = uniqueRowsCA(group1ca1id);
disp([num2str(size(uCA,1)),' control ca1 place cells'])
[uCA,~,~] = uniqueRowsCA(group2ca1id);
disp([num2str(size(uCA,1)),' pae ca1 place cells'])

[uCA,~,~] = uniqueRowsCA(group1ca3id);
disp([num2str(size(uCA,1)),' control ca3 place cells'])
[uCA,~,~] = uniqueRowsCA(group2ca3id);
disp([num2str(size(uCA,1)),' pae ca3 place cells'])


disp('...')
disp('control rats with ca1')
C=unique(extractBefore(unique(group1ca1id(:,1)),'_'));
for i=1:length(C)
    temp=uniqueRowsCA(group1ca1id);
    disp([C{i},' ',num2str(sum(contains(temp(:,1),C{i})))])
end

disp('...')
disp('lesion rats with ca1')
C=unique(extractBefore(unique(group2ca1id(:,1)),'_'));
for i=1:length(C)
    temp=uniqueRowsCA(group2ca1id);
    disp([C{i},' ',num2str(sum(contains(temp(:,1),C{i})))])
end

disp('...')
disp('control rats with ca3')
C=unique(extractBefore(unique(group1ca3id(:,1)),'_'));
for i=1:length(C)
    temp=uniqueRowsCA(group1ca3id);
    disp([C{i},' ',num2str(sum(contains(temp(:,1),C{i})))])
end

disp('...')
disp('lesion rats with ca3')
C=unique(extractBefore(unique(group2ca3id(:,1)),'_'));
for i=1:length(C)
    temp=uniqueRowsCA(group2ca3id);
    disp([C{i},' ',num2str(sum(contains(temp(:,1),C{i})))])
end



%% pop vec
% popvector(group1id,group1(:,end),'control');
% popvector(group2id,group2(:,end),'lesion');




% 
% visualize_cells(group1ca1id,'D:\Projects\HPCatn\figures\13_10_2019\controlca1')
% visualize_cells(group2ca1id,'D:\Projects\HPCatn\figures\13_10_2019\lesionca1')
% 
% visualize_cells(group1ca3id,'D:\Projects\HPCatn\figures\13_10_2019\controlca3')
% visualize_cells(group2ca3id,'D:\Projects\HPCatn\figures\13_10_2019\lesionca3')


group1ca1_trial_deviation = compute_trial_deviation(group1ca1id,group1ca1(:,contains(varnames,'runningdir')));

group2ca1_trial_deviation = compute_trial_deviation(group2ca1id,group1ca1(:,contains(varnames,'runningdir')));

group1ca3_trial_deviation = compute_trial_deviation(group1ca3id,group1ca3(:,contains(varnames,'runningdir')));

group2ca3_trial_deviation = compute_trial_deviation(group2ca3id,group2ca3(:,contains(varnames,'runningdir')));

group1ca1=[group1ca1,group1ca1_trial_deviation];
group2ca1=[group2ca1,group2ca1_trial_deviation];
group1ca3=[group1ca3,group1ca3_trial_deviation];
group2ca3=[group2ca3,group2ca3_trial_deviation];
% 
% varnames = [varnames,'trialdeviation'];



fig = figure;fig.Color = [1 1 1];
AllStats=stat_plot(group1ca3_trial_deviation,group2ca3_trial_deviation,{'Control','ATN Lesion'},{'Trial Deviation'},'plots',2)

fig = figure;fig.Color = [1 1 1];

AllStats=stat_plot(group1ca3(:,contains(varnames,'lap_perm_stability')),...
    group2ca3(:,contains(varnames,'lap_perm_stability')),{'Control','ATN Lesion'},{'Trial Stability'},'plots',2)



figure
AllStats=stat_plot(group1ca1(:,:,1),group2ca1(:,:,1),{'Control','ATN Lesion'},varnames)
figure
AllStats=stat_plot(group1ca3(:,:,1),group2ca3(:,:,1),{'Control','ATN Lesion'},varnames)

%% compile and save as csv
id=get_region_id([group1ca1id;group1ca3id;group2ca1id;group2ca3id],'F:\Projects\HPCatn\AnimalMetadata');

id(contains(id(:,end),'dg'),end) = cellstr(repmat('ca3',sum(contains(id(:,end),'dg')),1));

id = [id,[cellstr(repmat('control',size([group1ca1id;group1ca3id],1),1));...
    cellstr(repmat('lesion',size([group2ca1id;group2ca3id],1),1))]]

Rdata = [group1ca1;group1ca3;group2ca1;group2ca3];

Rdata = [num2cell(Rdata),id];

varnames = [varnames,{'session','tt','cell','area','group'}];
varnames = regexprep(varnames, '\W', '');
Rdata = cell2table(Rdata,'VariableNames',varnames);


writetable(Rdata,'F:\Projects\HPCatn\Rdata_hpcatn_sfn2019_control_lesion.csv')



%%
figure
plot(sort(group1ca3(:,55,1)),linspace(0,1,length(group1ca3(:,55,1))),'k')
hold on
plot(sort(group1ca3(:,55,2)),linspace(0,1,length(group1ca3(:,55,2))),'-.k')
plot(sort(group2ca3(:,55,1)),linspace(0,1,length(group2ca3(:,55,1))),'r')
plot(sort(group2ca3(:,55,2)),linspace(0,1,length(group2ca3(:,55,2))),'-.r')
xlabel(varnames{55})
box off
grid on

%% plot individual rats
clear colors
for i=1:length(control)
colors(:,i)=contains(group1ca1id(:,1),control{i})*i;
end
colors=sum(colors,2);

figure
plot(sort(group1ca1(colors==1,1,1)),linspace(0,1,length(group1ca1(colors==1,1,1))))
hold on
plot(sort(group1ca1(colors==2,1,1)),linspace(0,1,length(group1ca1(colors==2,1,1))))
plot(sort(group1ca1(colors==3,1,1)),linspace(0,1,length(group1ca1(colors==3,1,1))))

clear colors
for i=1:length(lesion)
colors(:,i)=contains(group2ca1id(:,1),lesion{i})*i;
end
colors=sum(colors,2);

plot(sort(group2ca1(colors==1,1,1)),linspace(0,1,length(group2ca1(colors==1,1,1))),'.-r')
plot(sort(group2ca1(colors==2,1,1)),linspace(0,1,length(group2ca1(colors==2,1,1))),'r')
xlabel('info content')
legend(control{1},control{2},control{3},lesion{1},lesion{2})

%% plot all measures to ppt
% for i=1:length(varnames)
%     AllStats=CDFplots(group1(:,i),group2(:,i),{'Control','ATN Lesion'},varnames(i),2);
%     toPPT(figure(1),'exportMode','matlab');
%     toPPT('setTitle',AllStats);
%     close all
% end

%% Locate Sessions for decoding
% clear data
% uniqueids=uniqueRowsCA(group1id);
% 
% uniquesessions=unique(group1id(:,1));
% for i=1:length(uniquesessions)
%     load(uniquesessions{i,1},'spikesID','linear_track')
%     
%     ensemble_rows=find(contains(uniqueids(:,1),uniquesessions{i,1}));
% 
%     pos=[linear_track.right{1,1}.dataspks(linear_track.right{1,1}.dataspks(:,6)==0,:);...
%         linear_track.left{1,1}.dataspks(linear_track.left{1,1}.dataspks(:,6)==0,:)];
%     [~,I]=sort(pos(:,1));
%     pos=pos(I,:);
%     
%     for k=1:length(ensemble_rows)
%         cells=find(contains(spikesID.TetrodeNum,uniqueids{ensemble_rows(k),2})...
%             & ismember(spikesID.CellNum,str2double(uniqueids(ensemble_rows(k),3))))';
%         spike_times{k,1}=sort([linear_track.right{1,cells}.dataspks(linear_track.right{1,cells}.dataspks(:,6)==1,1);...
%             linear_track.left{1,cells}.dataspks(linear_track.left{1,cells}.dataspks(:,6)==1,1)]);
%     end
%     pos_times=pos(:,1);
%     vel=pos(:,5);
%     head_angle=pos(:,4);
%     pos(:,[1,4:end])=[];
%     
%     
%     save(['D:\Projects\HPCatn\DecodingData\',uniquesessions{i,1}],'spike_times','pos','pos_times','vel','head_angle')
%     clear spike_times pos pos_times vel head_angle
% end


%%

function d = compute_trial_deviation(id,runningdir)
direction = {'left','right'};

parfor i=1:length(id)
    
    data=load(id{i,1},'linear_track','spikesID');
    
    idx = find_cells(data,str2double((extractBetween(id{i,2},'TT','.mat'))),str2double(id{i,3}));
    
    trials = vertcat(data.linear_track{1, 1}.(direction{runningdir(i)}){idx}.maps{:});
    
    d(i,1) = trial_deviation(trials);

end

end

function [group,groupid]=placefieldfilter(group,groupid,varnames)
% 1) Minimum peak firing rate of 2 Hz, 
% 2) Minimum field width of 8 cm, 
% 3) Maximum field width of 180 cm, 
% 4) at least 10 trials with consistent behavior. 
% 5) at least 100 spikes

groupid=groupid(group(:,contains(varnames,'PeakRate'))>=1 &...
    group(:,contains(varnames,'FieldWidth'))>=8 &...
    group(:,contains(varnames,'FieldWidth'))<=180 &...
    group(:,contains(varnames,'nlaps'))>=15 &...
    group(:,contains(varnames,'nSpikes'))>=100,:);

group=group(group(:,contains(varnames,'PeakRate'))>=1 &...
    group(:,contains(varnames,'FieldWidth'))>=8 &...
    group(:,contains(varnames,'FieldWidth'))<=180 &...
    group(:,contains(varnames,'nlaps'))>=15 &...
    group(:,contains(varnames,'nSpikes'))>=100,:,:);
end

function visualizecells(groupid,group)
cd('D:\Projects\HPCatn\ProcessedData')
for i=1:length(groupid)
    data=load(groupid{i,1});
%     close all
    postprocessFigures.main(data,{groupid{i,2},str2double(groupid(i,3))});
    set(gcf, 'Position', get(0, 'Screensize'));
%     pause(3)

    print(gcf,'-dpng', '-r90',...
        ['D:\Projects\HPCatn\figures\',group,...
        filesep,groupid{i,1},groupid{i,2},groupid{i,3},'.png'])
    close all
    
    
%     print(gcf,'-dpng', '-r90',...
%         ['G:\HPCatn\figures\',group,...
%         filesep,groupid{i,1},groupid{i,2},groupid{i,3},'.png'])
%     close all
%     


%     saveas(gcf,['D:\Projects\PAE_PlaceCell',filesep,'test.emf'])
%         p=postprocessFigures2(data,{groupid{i,2},str2double(groupid(i,3))});

% print('-dmeta',['D:\Projects\PAE_PlaceCell',filesep,'test.emf'])
end
end



function popvector(groupid,runningdir,group)
for i=1:length(groupid)
    data=load(groupid{i,1},'ratemap','spikesID');
    cells=find(contains(data.spikesID.TetrodeNum,groupid{i,2}) &...
        ismember(data.spikesID.CellNum,str2double(groupid{i,3})))';
    % rescale from 0 to 1 and save
    placecellmaps(i,:)=rescale(data.ratemap{cells,runningdir(i)},0,1);
    
    % save for opposite running direction
    if runningdir(i)==1
        oppositerunningmaps(i,:)=rescale(data.ratemap{cells,2},0,1);
    else
        oppositerunningmaps(i,:)=rescale(data.ratemap{cells,1},0,1);
    end


end

[~,I]=max(placecellmaps,[],2);
[~,I2]=sort(I);

placecellmaps=placecellmaps(I2,:);
oppositerunningmaps=oppositerunningmaps(I2,:);

fig=figure;fig.Color=[1 1 1];
subplot(1,2,1)
imagesc(placecellmaps)
ylabel('cells')
xlabel('position bins')
set(gca,'FontSize',12)

subplot(1,2,2)
imagesc(oppositerunningmaps)
ylabel('cells')
xlabel('position bins')
set(gca,'FontSize',12)


set(gcf, 'Position', get(0, 'Screensize'));
print(gcf,'-dpng', '-r400',...
    ['D:\Projects\HPCatn\figures\popvectors\',group,filesep,'popvec.png'])

savefig(['D:\Projects\HPCatn\figures\popvectors\',group,filesep,'popvec.fig'])
close all
end