% AnalyzeResults_HPCatn_inactivation
clear
data=compileResults('F:\Projects\HPCatn\ProcessedData');

rats_id = {'HPCatn06', 'HPCatn07'};

load('F:\Projects\HPCatn\successful_inactivation.mat')

for irats = 1:length(rats_id)
    load(['F:\Projects\HPCatn\AnimalMetadata\',rats_id{irats},'_metadata.mat'])
    cell_id = data.(rats_id{irats}).id;
    
    sessions=fieldnames(AnimalMetadata.RecordingLogs);
    for i=1:length(sessions)
        mazes{i}=AnimalMetadata.RecordingLogs.(sessions{i}).MazeTypes;
        notes{i}=AnimalMetadata.RecordingLogs.(sessions{i}).Notes;
    end
    
    % session types
    muscimol = sessions(contains(notes,'muscimol'));
    pbs = sessions(contains(notes,'pbs'));
    recovery = sessions(contains(notes,["6 hour"]));
    baseline = sessions(~contains(notes,["muscimol","pbs","6 hour","24 hour"]));
    
    fail_idx = contains(muscimol,extractBetween(failed_inactivation,'_','.mat'));
    pbs = [pbs;muscimol(fail_idx)];
    muscimol(fail_idx) = [];
    
    if irats==1
        idx(:,1)=contains(cell_id(:,1),sessions(contains(mazes,'track'))) & contains(cell_id(:,1),muscimol);
        idx(:,2)=contains(cell_id(:,1),sessions(contains(mazes,'track'))) & contains(cell_id(:,1),pbs);
        idx(:,3)=contains(cell_id(:,1),sessions(contains(mazes,'track'))) & contains(cell_id(:,1),recovery);
    elseif irats>1
        idx=[idx;[contains(cell_id(:,1),sessions(contains(mazes,'track'))) & contains(cell_id(:,1),muscimol),...
            contains(cell_id(:,1),sessions(contains(mazes,'track'))) & contains(cell_id(:,1),pbs),...
            contains(cell_id(:,1),sessions(contains(mazes,'track'))) & contains(cell_id(:,1),recovery)]];
    end
    
    clear notes mazes
end
varnames = [data.varnames,'runningdir','baseline','muscimol','pbs','recovery'];
measures = [data.(rats_id{1}).measures;data.(rats_id{2}).measures];
id = [data.(rats_id{1}).id;data.(rats_id{2}).id];



measures = measures(any(idx,2),:,:);
id = id(any(idx,2),:,:);
idx = idx(any(idx,2),:);

% unpack sessions data
r =size(measures,1);
measures = [[measures(:,:,1),ones(r,1)];...
    [measures(:,:,2),ones(r,1)+1];...
    [measures(:,:,3),ones(r,1)];...
    [measures(:,:,4),ones(r,1)+1]];

measures = [measures,[[ones(r,1);ones(r,1);zeros(r,1);zeros(r,1)],...
    [zeros(r,1);zeros(r,1);idx(:,1);idx(:,1)],...
    [zeros(r,1);zeros(r,1);idx(:,2);idx(:,2)],...
    [idx(:,3);idx(:,3);idx(:,3);idx(:,3)]]];

id = [id;id;id;id];

% locate brain areas
id=get_region_id(id,'F:\Projects\HPCatn\AnimalMetadata');

% filter for place cells
[measures_pc,id_pc]=placefieldfilter(measures,id,varnames);



% baseline vs ...
[ids_mus,current_measures_mus,comparison_data_mus]=baseline_vs('muscimol',measures,measures_pc,id,id_pc,varnames);
[ids_pbs,current_measures_pbs,comparison_data_pbs]=baseline_vs('pbs',measures,measures_pc,id,id_pc,varnames);


%% for R
varnames = [varnames,'condition','session','tt','cell','area'];

% baseline & muscimol
Rdata = [num2cell(current_measures_mus),cellstr(repmat('baseline_mus',size(current_measures_mus,1),1));...
    num2cell(comparison_data_mus),cellstr(repmat('muscimol',size(comparison_data_mus,1),1))];
Rdata = [Rdata,[ids_mus;ids_mus]];

% baseline & pbs
Rdata = [Rdata;[num2cell(current_measures_pbs),cellstr(repmat('baseline_pbs',size(current_measures_pbs,1),1));...
    num2cell(comparison_data_pbs),cellstr(repmat('pbs',size(comparison_data_pbs,1),1))],[ids_pbs;ids_pbs]];

%% add recovery data
idx = measures_pc(:,contains(varnames,'recovery'))==1;
ids = id_pc(idx,:);
current_measures = measures_pc(idx,:,:);

Rdata = [Rdata;[num2cell(current_measures),cellstr(repmat('recovery',size(current_measures,1),1)),ids]];

varnames = regexprep(varnames, '\W', '');

Rdata = cell2table(Rdata,'VariableNames',varnames);

Rdata.baseline(logical(Rdata.recovery))=0;

condition = Rdata.condition;

% code condition (baseline, drug, recovery)
Rdata.condition(logical(Rdata.baseline)) =...
    cellstr(repmat('baseline',sum(Rdata.baseline),1));
Rdata.condition(logical(Rdata.muscimol) | logical(Rdata.pbs)) =...
    cellstr(repmat('drug',sum(logical(Rdata.muscimol) | logical(Rdata.pbs)),1));
Rdata.condition(logical(Rdata.recovery)) =...
        cellstr(repmat('recovery',sum(Rdata.recovery),1));

% code drug (muscimol, pbs)    
Rdata.drug(logical(Rdata.muscimol) | strcmp(condition,'baseline_mus')) =...
    cellstr(repmat('muscimol',sum(Rdata.muscimol | strcmp(condition,'baseline_mus')),1));
Rdata.drug(logical(Rdata.pbs) | strcmp(condition,'baseline_pbs')) =...
    cellstr(repmat('pbs',sum(Rdata.pbs | strcmp(condition,'baseline_pbs')),1));

% date match to find recovery sessions
for i = 1:length(Rdata.session)
    cur_date = extractBetween(Rdata.session{i},'_S','.mat');
    all_date{i} = cur_date{1}(1:8);
end

mus_sessions = unique(Rdata.session(strcmp(Rdata.drug,'muscimol')));
for i = 1:length(mus_sessions)
    cur_date = extractBetween(mus_sessions{i},'_S','.mat');
    Rdata.drug(strcmp(all_date,cur_date{1}(1:8))) =...
        cellstr(repmat('muscimol',sum(strcmp(all_date,cur_date{1}(1:8))),1));
end

pbs_sessions = unique(Rdata.session(strcmp(Rdata.drug,'pbs')));
for i = 1:length(pbs_sessions)
    cur_date = extractBetween(pbs_sessions{i},'_S','.mat');
    Rdata.drug(strcmp(all_date,cur_date{1}(1:8))) =...
        cellstr(repmat('pbs',sum(strcmp(all_date,cur_date{1}(1:8))),1));
end

Rdata.drug(ismember(Rdata.session,unique(Rdata.session(strcmp(Rdata.drug,'muscimol'))))) = ...
    cellstr(repmat('muscimol',sum(ismember(Rdata.session,unique(Rdata.session(strcmp(Rdata.drug,'muscimol'))))),1));

Rdata.drug(ismember(Rdata.session,unique(Rdata.session(strcmp(Rdata.drug,'pbs'))))) = ...
    cellstr(repmat('pbs',sum(ismember(Rdata.session,unique(Rdata.session(strcmp(Rdata.drug,'pbs'))))),1));

Rdata(contains(Rdata.session,'HPCatn06_S20190301165853.mat'),:) = [];

% writetable(Rdata,'F:\Projects\HPCatn\Rdata_hpcatn_sfn2019_inactivation.csv')

idx = contains(Rdata.condition,'baseline') & contains(Rdata.drug,'pbs');
uniqueRowsCA([Rdata.session(idx),Rdata.tt(idx),Rdata.cell(idx)])

idx = contains(Rdata.condition,'baseline') & contains(Rdata.drug,'muscimol');
uniqueRowsCA([Rdata.session(idx),Rdata.tt(idx),Rdata.cell(idx)])

idx = contains(Rdata.condition,'recovery') & contains(Rdata.drug,'pbs');
uniqueRowsCA([Rdata.session(idx),Rdata.tt(idx),Rdata.cell(idx)])

idx = contains(Rdata.condition,'recovery') & contains(Rdata.drug,'muscimol');
uniqueRowsCA([Rdata.session(idx),Rdata.tt(idx),Rdata.cell(idx)])

%%
for i = 1:length(Rdata.session)
   load(Rdata.session{i},'binned_vel');
   if contains(Rdata.condition{i},'baseline')
       Rdata.vel(i) = nanmean(binned_vel{1});
   elseif contains(Rdata.condition{i},'drug')
       Rdata.vel(i) = nanmean(binned_vel{2});
   elseif contains(Rdata.condition{i},'recovery')
       Rdata.vel(i) = nanmean(binned_vel{1});
   end
end
%%
% tempid = [Rdata.session(contains(Rdata.condition,'drug') & contains(Rdata.drug,'muscimol')),...
%     Rdata.tt(contains(Rdata.condition,'drug') & contains(Rdata.drug,'muscimol')),...
%     Rdata.cell(contains(Rdata.condition,'drug') & contains(Rdata.drug,'muscimol'))];
% 
% visualize_cells(uniqueRowsCA(tempid),'F:\Projects\HPCatn\figures\place_cells\muscimol')
% 
% 
% tempid = [Rdata.session(contains(Rdata.condition,'drug') & contains(Rdata.drug,'pbs')),...
%     Rdata.tt(contains(Rdata.condition,'drug') & contains(Rdata.drug,'pbs')),...
%     Rdata.cell(contains(Rdata.condition,'drug') & contains(Rdata.drug,'pbs'))];
% 
% visualize_cells(uniqueRowsCA(tempid),'F:\Projects\HPCatn\figures\place_cells\pbs')
%%
var = 1;
for i = 1:length(varnames)
    figure
    AllStats{i,1}=stat_plot(Rdata.(varnames{i})(contains(Rdata.condition,'drug') & contains(Rdata.drug,'pbs')),...
        Rdata.(varnames{i})(contains(Rdata.condition,'drug') & contains(Rdata.drug,'muscimol')),...
        {'pbs','muscimol'},varnames{i})
end

for i = 1:length(varnames)
    figure
    AllStats{i,1}=stat_plot(Rdata.(varnames{i})(contains(Rdata.condition,'recovery') & contains(Rdata.drug,'pbs')),...
        Rdata.(varnames{i})(contains(Rdata.condition,'recovery') & contains(Rdata.drug,'muscimol')),...
        {'pbs_recovery','muscimol_recovery'},varnames{i})
end

% drug session / pbs vs muscimol
AllStats=stat_plot(Rdata.(varnames{var})(contains(Rdata.condition,'drug') & contains(Rdata.drug,'pbs')),...
    Rdata.(varnames{var})(contains(Rdata.condition,'drug') & contains(Rdata.drug,'muscimol')),...
    {'pbs','muscimol'},varnames{var})

% baseline session / pbs vs muscimol
AllStats=stat_plot(Rdata.(varnames{var})(contains(Rdata.condition,'baseline') & contains(Rdata.drug,'pbs')),...
    Rdata.(varnames{var})(contains(Rdata.condition,'baseline') & contains(Rdata.drug,'muscimol')),...
    {'pbs','muscimol'},varnames{var},'plottype','beeswarm')

% recovery session / pbs vs muscimol
AllStats=stat_plot(Rdata.(varnames{var})(contains(Rdata.condition,'recovery') & contains(Rdata.drug,'pbs')),...
    Rdata.(varnames{var})(contains(Rdata.condition,'recovery') & contains(Rdata.drug,'muscimol')),...
    {'pbs','muscimol'},varnames{var},'plottype','beeswarm')

% baseline vs recovery session / pbs vs pbs
AllStats=stat_plot(Rdata.(varnames{var})(contains(Rdata.condition,'baseline') & contains(Rdata.drug,'pbs')),...
    Rdata.(varnames{var})(contains(Rdata.condition,'recovery') & contains(Rdata.drug,'pbs')),...
    {'pbs_baseline','pbs recovery'},varnames{var})

% baseline vs recovery session / muscimol vs muscimol
AllStats=stat_plot(Rdata.(varnames{var})(contains(Rdata.condition,'baseline') & contains(Rdata.drug,'muscimol')),...
    Rdata.(varnames{var})(contains(Rdata.condition,'recovery') & contains(Rdata.drug,'muscimol')),...
    {'pbs_baseline','pbs recovery'},varnames{var})

% baseline vs drug session / muscimol vs muscimol
AllStats=stat_plot(Rdata.(varnames{var})(contains(Rdata.condition,'baseline') & contains(Rdata.drug,'muscimol')),...
    Rdata.(varnames{var})(contains(Rdata.condition,'drug') & contains(Rdata.drug,'muscimol')),...
    {'pbs_baseline','pbs recovery'},varnames{var})


% %% look at session xy to see if rat was affected by drug
% session_to_test = unique(Rdata.session(strcmp(Rdata.drug,'muscimol') & strcmp(Rdata.condition,'drug')));
% clear answer_
% for i = 1:length(session_to_test)
%     load(session_to_test{i},'linear_track')
%     if size(linear_track,2)~=2
%         answer_{i} = 'r';
%         continue
%     end
%     
%     fig = figure; fig.Color = [1 1 1];
%     subplot(2,1,1)
%     plot(linear_track{1, 1}.nonlinearFrames(:,2),linear_track{1, 1}.nonlinearFrames(:,3),'.k')
%     axis image
%     grid on
%     title(session_to_test{i})
% 
%     subplot(2,1,2)
%     plot(linear_track{1, 2}.nonlinearFrames(:,2),linear_track{1, 2}.nonlinearFrames(:,3),'.k')
%     axis image
%     grid on
%     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.4958 0.0314 0.5083 0.9686]);
%     answer_{i} = input('Successful inactivation? (''y'' ''n'' ''r''): ','s');
%     close 
% end
% successful_inactivation = session_to_test(strcmp(answer_,'y'))
% 
% failed_inactivation = session_to_test(strcmp(answer_,'n'))

% table2array(Rdata); 
% Rdata()


% var = 1;
% fig=figure; fig.Color = [1 1 1];
% for i=1:size(current_measures,1)
%     if current_measures(i,var)>muscimol_comparison(i,var) 
%         colors = [0.470000 0.810000 0.940000];
%     else
%         colors = [1.000000 0.700000 0.280000];
%     end
% plot([1,2],[current_measures(i,var),muscimol_comparison(i,var)],'Color',colors);hold on
% end
% scatter(ones(size(current_measures,1),1),current_measures(:,var),'k','Filled')
% scatter(ones(size(current_measures,1),1)+1,muscimol_comparison(:,var),'r','Filled')
% xlim([.5 2.5])
% grid on
% box off
% ylabel(varnames{var})
% set(gca,'XTick',[1 2]);
% set(gca,'XTickLabel',{'Baseline',comparison})
% set(gca,'FontSize',10,'FontWeight','bold','LineWidth',1,'box','off')

figure;
distributionPlot([current_measures(:,var),comparison_data(:,var)],...
    'showMM',3,'xNames',{'Baseline',comparison},'histOpt',1);hold on


%          subplot(3,4,12),distributionPlot(chi2rnd(3,1000,1),'histOri','right','color','r','widthDiv',[2 2],'showMM',0)
%          subplot(3,4,12),distributionPlot(chi2rnd(5,1000,1),'histOri','left','color','b','widthDiv',[2 1],'showMM',0)
%          
AllStats=stat_plot(current_measures(:,var),comparison_data(:,var),{'Baseline',comparison},varnames{var},'plottype','beeswarm')

% distributionPlot([current_measures(:,var),muscimol_comparison(:,var)],'addSpread',true,'showMM',3) 

AllStats=stat_plot(measures_pc(measures_pc(:,contains(varnames,'baseline'))==1 & measures_pc(:,contains(varnames,'recovery'))==0,:,:),...
    comparison_data,{'baseline','ATN Inactivation'},varnames)


figure
AllStats=stat_plot(measures(measures(:,contains(varnames,'baseline'))==1,:),...
    measures(measures(:,contains(varnames,'muscimol'))==1,:),{'baseline','ATN Inactivation'},varnames)

figure
AllStats=stat_plot(measures(measures(:,contains(varnames,'baseline'))==1,:),...
    measures(measures(:,contains(varnames,'pbs'))==1,:),{'baseline','ATN pbs'},varnames)

figure
AllStats=stat_plot(measures(measures(:,contains(varnames,'recovery'))==1,:),...
    measures(measures(:,contains(varnames,'muscimol'))==1,:),{'recovery','ATN Inactivation'},varnames)

figure
AllStats=stat_plot(measures(measures(:,contains(varnames,'recovery'))==1,:),...
    measures(measures(:,contains(varnames,'pbs'))==1,:),{'recovery','ATN pbs'},varnames)

figure
AllStats=stat_plot(measures(measures(:,contains(varnames,'baseline'))==1,:),...
    measures(measures(:,contains(varnames,'recovery'))==1,:),{'baseline','ATN recovery'},varnames)








function [ids,current_measures,comparison_data]=baseline_vs(comparison,measures,measures_pc,id,id_pc,varnames)
% baseline vs ...
idx = measures_pc(:,contains(varnames,'baseline'))==1 & measures_pc(:,contains(varnames,'recovery'))==0;
ids = id_pc(idx,:);

current_measures = measures_pc(idx,:,:);

runningdir = measures_pc(measures_pc(:,contains(varnames,'baseline'))==1,...
    contains(varnames,'runningdir'));

for i=1:length(ids)
    idx = strcmp(id(:,1),ids(i,1)) &...
        strcmp(id(:,2),ids(i,2)) &...
        strcmp(id(:,3),ids(i,3)) &...
        strcmp(id(:,4),ids(i,4)) &...
        measures(:,contains(varnames,comparison))==1 &...
        measures(:,contains(varnames,'runningdir'))==runningdir(i);
    if ~any(idx)
        comparison_data(i,:) = NaN(1,size(measures,2));
        cannot_locate(i,:)=1;
    else
        comparison_data(i,:) = measures(idx,:);
        cannot_locate(i,:)=0;
    end
end
ids(logical(cannot_locate),:)=[];
current_measures(logical(cannot_locate),:)=[];
comparison_data(logical(cannot_locate),:)=[];
end



%%
function [group,groupid]=placefieldfilter(group,groupid,varnames)
% 1) Minimum peak firing rate of 2 Hz, 
% 2) Minimum field width of 8 cm, 
% 3) Maximum field width of 180 cm, 
% 4) at least 15 trials with consistent behavior. 
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