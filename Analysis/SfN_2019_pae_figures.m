% SfN_2019_pae_figures

%cell examples

% control 
% ca1
linear_track_example('LEM3116_S20180720112057.mat',2,3,'right')
linear_track_example('LEM3216_S20190718101925.mat',15,5,'right')
linear_track_example('LEM3216_S20190718101925.mat',3,1,'left')
linear_track_example('LEM3216_S20190718101925.mat',6,2,'left')


cylinder_example('LEM3216_S20190710184004.mat',2,6,2)
cylinder_example('LEM3216_S20190710184004.mat',2,4,3)
cylinder_example('LEM3216_S20190712114116.mat',6,2,2)
cylinder_example('LEM3216_S20190712114116.mat',6,2,3)
cylinder_example('LEM3216_S20190712114116.mat',6,7,2)
cylinder_example('LEM3216_S20190712114116.mat',6,16,3)
cylinder_example('LEM3216_S20190717174058.mat',16,10,2)


% ca3
linear_track_example('LEM3116_S20180806104145.mat',8,1,'left') 
linear_track_example('LEM3116_S20180807105247.mat',8,4,'left') 
linear_track_example('LEM3216_S20190807105102.mat',2,12,'right') 


cylinder_example('LEM3116_S20180806104145.mat',2,1,3)
cylinder_example('LEM3116_S20180806104145.mat',8,1,3)
cylinder_example('LEM3116_S20180807105247.mat',2,13,3)
cylinder_example('LEM3116_S20180807105247.mat',2,12,2)
cylinder_example('LEM3116_S20180807105247.mat',3,2,3)
cylinder_example('LEM3116_S20180807105247.mat',4,6,2)
cylinder_example('LEM3116_S20180807105247.mat',6,2,3)
cylinder_example('LEM3216_S20190807105102.mat',2,16,2) 



% dg
linear_track_example('LEM3116_S20180727123953.mat',5,2,'left')
linear_track_example('LEM3116_S20180731114926.mat',5,1,'left')
linear_track_example('LEM3216_S20190726184722.mat',13,118,'right')


cylinder_example('LEM3116_S20180727123953.mat',5,2,2)
cylinder_example('LEM3116_S20180727123953.mat',5,2,3)
cylinder_example('LEM3116_S20180802100324.mat',5,6,2)
% cylinder_example('LEM3116_S20180803103321.mat',2,1,3)
cylinder_example('LEM3116_S20180815132530.mat',3,2,3)





%% sizing
set(findall(gcf,'-property','FontSize'),'FontSize',27)
set(findall(gcf,'-property','FontSize'),'FontSize',16)

%% pae
% ca1
linear_track_example("LS19_S20170508205640.mat",4,5,'right')
linear_track_example("LS19_S20170523165204.mat",6,2,'left')
linear_track_example("LEM3124_S20190302163218.mat",4,7,'left')


cylinder_example('LS19_S20170523215606.mat',1,3,2)
cylinder_example('LEM3124_S20190228180140.mat',4,21,2)
cylinder_example('LEM3124_S20190228180140.mat',11,7,2)
cylinder_example('LEM3124_S20190228180140.mat',11,9,3)
cylinder_example('LEM3124_S20190301133246.mat',4,14,2)

% ca3
linear_track_example("LEM3124_S20190307163732.mat",13,4,'left')
linear_track_example("LEM3124_S20190311161641.mat",2,9,'left')
linear_track_example("LEM3124_S20190315160017.mat",3,8,'right')

cylinder_example('LEM3124_S20190307163732.mat',4,1,2)
cylinder_example('LEM3124_S20190307163732.mat',7,4,2)
cylinder_example('LEM3124_S20190307163732.mat',7,6,2)


% dg
linear_track_example('LEM3246_S20190629164240.mat',8,6,'right')
linear_track_example('LEM3246_S20190701153855.mat',1,11,'right')
linear_track_example('LEM3124_S20190307132509.mat',11,5,'left')

cylinder_example('LEM3124_S20190307132509.mat',7,10,2)
cylinder_example('LEM3124_S20190307132509.mat',7,4,2)
cylinder_example('LEM3124_S20190307132509.mat',7,13,2)


% 
% linear_track_example('LEM3216_S20190718101925.mat',8,11,'right')
% set(findall(gcf,'-property','FontSize'),'FontSize',13)
% 
% linear_track_example('LEM3216_S20190718101925.mat',6,2,'right')
% set(findall(gcf,'-property','FontSize'),'FontSize',13)
% 
% linear_track_example('LEM3216_S20190718101925.mat',3,1,'right')
% set(findall(gcf,'-property','FontSize'),'FontSize',13)
% 
% 
% linear_track_example('LEM3216_S20190718101925.mat',3,1,'right')
% set(findall(gcf,'-property','FontSize'),'FontSize',13)


%% Theta AutoCorr
 
sessions = uniqueRowsCA([track_cylinder.session,track_cylinder.tt,track_cylinder.cell,track_cylinder.group]);
clear autocorr
% for i = 1:length(sessions)
%     disp([sessions{i,:}])
%     data = load(sessions{i,1},'spikesID','thetaautocorr');
%     idx = find_cells(data,str2double(extractBetween(sessions{i,2},'TT','.mat')),str2double(sessions{i,3}));
%     
%     if contains(track_cylinder.mazetype,'cylinder')
%         autocorr(i,:) = data.thetaautocorr{idx,3};
%     else
%         autocorr(i,:) = data.thetaautocorr{idx,track_cylinder.runningdir(i)};
%     end
% end
for i = 1:length(sessions)
    disp([sessions{i,:}])
    data = load(sessions{i,1},'spikesID','Spikes');
    idx = find_cells(data,str2double(extractBetween(sessions{i,2},'TT','.mat')),str2double(sessions{i,3}));
    
    [thetaindex(i,1),thetapeak(i,1),cor(i,:),~]=thetamodulation(data.Spikes{idx});


end
autocorr = cor;

fig = figure;fig.Color = [1 1 1];
[~,idx]=sort(thetaindex(contains(sessions(:,end),'control'),:),'descend');
control_autocorr = autocorr(contains(sessions(:,end),'control'),:);
control_autocorr = control_autocorr(idx,:);
imagesc(control_autocorr);
colormap(viridis)
axis xy
hold on
plot(1:size(control_autocorr,2),rescale(nanmean(control_autocorr),1,size(control_autocorr,1)),'w','linewidth',2);
set(gca,'XMinorTick','on','YMinorTick','off',...
    'LineWidth',1,'box','off','FontSize',20,'FontWeight','bold','XTickLabels',[])

fig = figure;fig.Color = [1 1 1];
[~,idx]=sort(thetaindex(contains(sessions(:,end),'pae'),:),'descend');
control_autocorr = autocorr(contains(sessions(:,end),'pae'),:);
control_autocorr = control_autocorr(idx,:);
imagesc(control_autocorr);
colormap(viridis)
axis xy
hold on
plot(1:size(control_autocorr,2),rescale(nanmean(control_autocorr),1,size(control_autocorr,1)),'w','linewidth',2);
set(gca,'XMinorTick','on','YMinorTick','off',...
    'LineWidth',1,'box','off','FontSize',20,'FontWeight','bold','XTickLabels',[])

AllStatsca1=stat_plot(thetaindex(sessions(:,end) == 'control'),thetaindex(sessions(:,end) == 'pae'),{'Sacc','PAE'},'thetaindex')

% fig = figure;fig.Color = [1 1 1];
% [~,idx]=sort(track_cylinder.thetaindex(contains(track_cylinder.group,'control'),:),'descend');
% control_autocorr = autocorr(contains(track_cylinder.group,'control'),:);
% control_autocorr = control_autocorr(idx,:);
% imagesc(control_autocorr);
% colormap(viridis)
% axis xy
% hold on
% plot(1:size(control_autocorr,2),rescale(nanmean(control_autocorr),1,size(control_autocorr,1)),'w','linewidth',2);
% set(gca,'XMinorTick','on','YMinorTick','off',...
%     'LineWidth',1,'box','off','FontSize',20,'FontWeight','bold','XTickLabels',[])
% 
% fig = figure;fig.Color = [1 1 1];
% [~,idx]=sort(track_cylinder.thetaindex(contains(track_cylinder.group,'pae'),:),'descend');
% pae_autocorr = autocorr(contains(track_cylinder.group,'pae'),:);
% pae_autocorr = pae_autocorr(idx,:);
% imagesc(pae_autocorr);
% colormap(viridis)
% axis xy
% hold on
% plot(1:size(pae_autocorr,2),rescale(nanmean(pae_autocorr),1,size(pae_autocorr,1)),'w','linewidth',2);
% set(gca,'XMinorTick','on','YMinorTick','off',...
%     'LineWidth',1,'box','off','FontSize',20,'FontWeight','bold','XTickLabels',[])

AllStatsca1=stat_plot(track_cylinder.thetaindex(track_cylinder.group == 'control'),track_cylinder.thetaindex(track_cylinder.group == 'pae'),{'Sacc','PAE'},'thetaindex')

%%
function cylinder_example(session,tt,cell,ns)
data = load(session);

[i]=find_cells(data,tt,cell);

tetrode=strsplit(data.spikesID.paths{i},filesep);
tetrode=tetrode{end};
trodeID=str2double(extractBetween(tetrode,'TT','.'));
colorcode = 'r';

fig=figure('Name',[data.rat,'  ',data.sessionID,'  ',tetrode,' Cell: ',num2str(data.spikesID.CellNum(i))],'NumberTitle','off');
fig.Color = [1 1 1];
p = panel(fig);

p.pack(4, 1);
% set margins
p.de.margin = 6;

p(1, 1).select();
ax = gca;
[data_video_spk,data_video_nospk]=createframes_w_spikebinary(data,ns,i);
postprocessFigures.spikesonpath_2d(ax,data_video_spk,data.lfp.ts,data.lfp.theta_phase(trodeID,:),colorcode)
title('')
disp(sprintf('IC: %4.2f  %4.2f hz',...
    data.measures(contains(data.varnames,["InformationContent","PeakRate"]))))


p(2, 1).select();
ax = gca;
postprocessFigures.ratemaps_2d(ax,data.ratemap{i,ns+1},data.measures(i,:,ns+1),data.varnames)
title('')


binside=mean([range(data_video_nospk(:,2))/length(data.ratemap{i,ns+1}),...
    range(data_video_nospk(:,3))/length(data.ratemap{i,ns+1})]);

results=pass_index(data_video_nospk(:,1),data_video_nospk(:,2:3),...
    data_video_spk(data_video_spk(:,6)==1,1),...
    [data.lfp.ts(data.lfp.ts>=data.events(1,ns) & data.lfp.ts<=data.events(2,ns))]',...
    [data.lfp.signal(trodeID,data.lfp.ts>=data.events(1,ns) & data.lfp.ts<=data.events(2,ns))]',...
    'plots',0,'method','place','binside',round(binside),'sample_along','arc_length');

            
bins=length(data.ratemap{i,ns+1});
xedge=linspace(-1,1,bins+1);
phaseedge=linspace(0,720,bins*4);


phase_spk = interp1(data.lfp.ts,data.lfp.theta_phase(trodeID,:),data_video_spk(data_video_spk(:,6)==1,1));
x_spk = interp1(results.ts,results.pass_index,data_video_spk(data_video_spk(:,6)==1,1));
spkmap=histcounts2([x_spk;x_spk],[phase_spk;phase_spk+2*pi]*180/pi,xedge,phaseedge);


phase_occ = interp1(data.lfp.ts,data.lfp.theta_phase(trodeID,:),data_video_spk(data_video_spk(:,6)==0,1));
x_occ = interp1(results.ts,results.pass_index,data_video_spk(data_video_spk(:,6)==0,1));
occ=histcounts2([x_occ;x_occ],[phase_occ;phase_occ+2*pi]*180/pi,xedge,phaseedge);
occ = occ/data.samplerate;

phasemap=spkmap./occ;

phasemap(isnan(phasemap)) = 0;
phasemap(isinf(phasemap)) = 0;

h=4;
phase_size = size(phasemap,2);
myfilter = fspecial('gaussian',[4 24]*h, h);
phasemap = imfilter([phasemap,phasemap,phasemap],myfilter,'replicate');
phasemap = phasemap(:,phase_size:phase_size*2);

p(3, 1).select();
scatter([x_spk;x_spk],[phase_spk;phase_spk+2*pi],20,'Filled','k');
axis tight
axis square
axis off

p(4, 1).select();
pcolor(flipud(rot90(phasemap)));shading flat;box off;axis off;axis tight;
colormap(viridis(255))
axis square
axis on
ylabel('\theta Phase')
xlabel('')
% set(gca,'XTick',linspace(ax.XLim(1),ax.XLim(2),3),'XTickLabels',[-1,0,1])
set(gca,'XTickLabels',[])
ax = gca;
set(gca,'YTick',linspace(ax.YLim(1),ax.YLim(2),3),'YTickLabels',{'0','2\pi','4\pi'})
            
set(gcf,'OuterPosition',[1050 34 369 1017])

end

function linear_track_example(session,tt,cell,track_direction)
data = load(session);

[i]=find_cells(data,tt,cell);

tetrode=strsplit(data.spikesID.paths{i},filesep);
tetrode=tetrode{end};
trodeID=str2double(extractBetween(tetrode,'TT','.'));
colorcode = 'r';


fig=figure('Name',[data.rat,'  ',data.sessionID,'  ',tetrode,' Cell: ',num2str(data.spikesID.CellNum(i))],'NumberTitle','off');
fig.Color = [1 1 1];
p = panel(fig);


p.pack(4, 1);
% set margins
p.de.margin = 3;

p(1, 1).select();
ax = gca;
ns = 1;
postprocessFigures.unlinearpath(ax,data.linear_track{round(ns/2)}.nonlinearFrames,...
    data.linear_track{round(ns/2)}.(track_direction){1,i}.laps,...
    data.linear_track{round(ns/2)}.(track_direction){1,i}.dataspks,data.mazetypes,...
    data.lfp.ts,data.lfp.theta_phase(trodeID,:),colorcode)
title('')

p(2, 1).select();
ax = gca;
if contains(track_direction,'left')
    ns = 2;
end
ratemaps_(ax,data.linear_track{round(ns/2)}.(track_direction){1,i}.maps,data.ratemap{i,ns},...
    data.measures(i,:,ns),data.varnames)
set(gca,'XTick',[])
ylabel('laps')
xlabel('')
ax = gca;
set(gca,'YTick',linspace(ax.YLim(1),ax.YLim(2),4),'YTickLabels',floor(linspace(ax.YLim(1),ax.YLim(2),4)))


p(3, 1).select();
ax = gca;
phase_by_pos_(ax,data.linear_track{round(ns/2)}.(track_direction){1,i}.dataspks,...
    data.lfp.ts,data.lfp.theta_phase(trodeID,:),...
    data.ratemap{i,ns},data.linear_track{round(ns/2)}.left{i})
title('')

% PHASE MAP
p(4, 1).select();
ax = gca;
postprocessFigures.phase_map(ax,data.linear_track{round(ns/2)}.(track_direction){1,i}.dataspks(data.linear_track{round(ns/2)}.(track_direction){1,i}.dataspks(:,6)==0,1),...
    data.linear_track{round(ns/2)}.(track_direction){1,i}.dataspks(data.linear_track{round(ns/2)}.(track_direction){1,i}.dataspks(:,6)==0,2),...
    data.linear_track{round(ns/2)}.(track_direction){1,i}.dataspks(data.linear_track{round(ns/2)}.(track_direction){1,i}.dataspks(:,6)==1,1),...
    data.lfp,data.ratemap{i,ns},trodeID,data.samplerate)
axis on
ylabel('\theta Phase')
xlabel('Distance')
set(gca,'YTick',[1 18.5 36],'YTickLabels',{'0','2\pi','4\pi'})

set(gca,'XTick',linspace(1,40,4),'XTickLabels',{'0','40','80','120'})


set(gcf,'OuterPosition',[1088 226 509 714])

end

function ratemaps_(ax,lapmaps,ratemap,measures,varnames)
% PLOT EACH LAPS RATEMAP WITH THE OVERALL FIRING RATE
% SUPERIMPOSED
if ~iscell(lapmaps)
    return
end
laps=reshape([lapmaps{:}],[],length(lapmaps));
imagesc(flipud(imrotate(laps,90)));hold on;set(gca,'Ydir','Normal')

% plot(rescale(ratemap,1,size(laps,2)),'LineWidth',2, 'color','w');
plot(rescale(nanmean(flipud(imrotate(laps,90))),1,size(laps,2)),'LineWidth',2, 'color','w');

set(gca,'XTick',linspace(1,120,7),'XTickLabel',0:360/6:360,...
    'FontSize',12)
xlabel('Distance (cm)')
axis tight
hold on;box off;
colormap(ax,viridis(255))
title(sprintf('IC: %4.2f  %4.2f hz',...
    measures(contains(varnames,["InformationContent","PeakRate"]))))
end

function phase_by_pos_(ax,dataspks,lfp_ts,theta_phase,ratemap,trackinfo)
if isempty(dataspks)
    return
end

% PHASE BY POSITION
x=rescale(dataspks(:,2),0,1);

phase=interp1(lfp_ts,theta_phase,dataspks(logical(dataspks(:,6)),1)','linear');
scatter([x(logical(dataspks(:,6)));x(logical(dataspks(:,6)))],[phase';phase'+2*pi]*180/pi,15,'Filled','k');

box off;axis off
xlim([min(x) max(x)])
set(ax,'YTick',0:360:720,'YTickLabel',[{'0'},{'2\pi'},{'4\pi'}])

if isfield(trackinfo,'fields')
    for f=1:length(trackinfo.fields)
        rho(f)=trackinfo.fields{f}.ThPrecess.circLinCorr;
        pval(f)=trackinfo.fields{f}.ThPrecess.pval;
    end
    a3=[rho;pval];
    title(sprintf(repmat(' rho:%4.2f p:%4.2f,',1,length(rho)),a3(:)'))
    clear rho pval
end
end