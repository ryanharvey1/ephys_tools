
addpath('F:\ClarkP30_Recordings\AnimalMetaData')
rats=dir(fullfile('F:\ClarkP30_Recordings\AnimalMetaData','*.mat'));
rats={rats.name};

for i=1:length(rats)
    load(rats{i})
    sess=fieldnames(AnimalMetadata.RecordingLogs);
    for s=1:length(sess)
        idx=contains(groupid(:,1),[extractBefore(rats{i},'_'),'_',sess{s}]);
        if sum(idx) == 0
            continue
        end
        currentsess=groupid(idx,:);
        if AnimalMetadata.RecordingLogs.(sess{s}).MazeTypes == "Cylinder"
            groupid(idx,4)={'baseline'};
        elseif AnimalMetadata.RecordingLogs.(sess{s}).MazeTypes == "Cylinder,Cylinder"
            groupid(idx,4)={'stability'};
        elseif AnimalMetadata.RecordingLogs.(sess{s}).MazeTypes == "Cylinder,Cylinder,Cylinder"
            groupid(idx,4)={'dim_cue_rotation'};
        elseif AnimalMetadata.RecordingLogs.(sess{s}).MazeTypes == "Cylinder,Cylinder,Cylinder,Cylinder" || ...
                AnimalMetadata.RecordingLogs.(sess{s}).MazeTypes == "Cylinder,Cylinder,Cylinder,Cylinder,box" || ...
                AnimalMetadata.RecordingLogs.(sess{s}).MazeTypes == "Cylinder,Cylinder,Cylinder,Cylinder,circ track" || ...
                AnimalMetadata.RecordingLogs.(sess{s}).MazeTypes == "Cylinder,Cylinder,Cylinder,Cylinder,Cylinder" 
            groupid(idx,4)={'stability_rotation'};
        elseif AnimalMetadata.RecordingLogs.(sess{s}).MazeTypes == "Cylinder,box,Cylinder"
            groupid(idx,4)={'delta_context'};
        end
    end
end

ang_diff = zeros(length(groupid),4);
for i = 1:length(groupid)
    
    if isempty(groupid(i,4))
        
        ang_diff(i,1:5) = nan;
    elseif strcmp(groupid(i,4),'baseline')
        
        ang_diff(i,1:5) = nan;
        
    elseif strcmp(groupid(i,4),'stability')
        
        PD_s1 = deg2rad(dir_measures_all(i,contains(vars,'preferred_Direction'),1)); 
        PD_s2 = deg2rad(dir_measures_all(i,contains(vars,'preferred_Direction'),2)); 
        ang_diff(i,1) = nan;
        ang_diff(i,2) = wrapTo360(rad2deg(circ_dist(PD_s1,PD_s2)));
        ang_diff(i,3) = nan; 
        ang_diff(i,4) = nan; 
        ang_diff(i,5) = nan;
        
    elseif strcmp(groupid(i,4),'dim_cue_rotation') || strcmp(groupid(i,4),'delta_context')
        
        PD_s1 = deg2rad(dir_measures_all(i,contains(vars,'preferred_Direction'),1)); 
        PD_s2 = deg2rad(dir_measures_all(i,contains(vars,'preferred_Direction'),2)); 
        PD_s3 = deg2rad(dir_measures_all(i,contains(vars,'preferred_Direction'),3)); 
        ang_diff(i,1) = nan;
        ang_diff(i,2) = nan; % rotation
        ang_diff(i,3) = wrapTo360(rad2deg(circ_dist(PD_s1,PD_s2))); % rotation
        ang_diff(i,4) = wrapTo360(rad2deg(circ_dist(PD_s2,PD_s3))); % return
        ang_diff(i,5) = wrapTo360(rad2deg(circ_dist(PD_s1,PD_s3))); % baseline versus return

        
    elseif strcmp(groupid(i,4),'stability_rotation')
        PD_s1 = deg2rad(dir_measures_all(i,contains(vars,'preferred_Direction'),1)); 
        PD_s2 = deg2rad(dir_measures_all(i,contains(vars,'preferred_Direction'),2)); 
        PD_s3 = deg2rad(dir_measures_all(i,contains(vars,'preferred_Direction'),3));
        PD_s4 = deg2rad(dir_measures_all(i,contains(vars,'preferred_Direction'),4));
        ang_diff(i,1) = nan; 
        ang_diff(i,2) = wrapTo360(rad2deg(circ_dist(PD_s1,PD_s2))); % stability
        ang_diff(i,3) = wrapTo360(rad2deg(circ_dist(PD_s2,PD_s3))); % rotation
        ang_diff(i,4) = wrapTo360(rad2deg(circ_dist(PD_s3,PD_s4))); % return
        ang_diff(i,5) = wrapTo360(rad2deg(circ_dist(PD_s1,PD_s4))); % baseline versus return
    end
end


WT_ang_diff = ang_diff(HD_idx_final & ~genotype,2);
WT_ang_diff(isnan(WT_ang_diff))=[];
wrapTo360(rad2deg(circ_mean(deg2rad(WT_ang_diff))))
circ_r(deg2rad(WT_ang_diff))
circ_kappa(deg2rad(WT_ang_diff))

Tg_ang_diff = ang_diff(HD_idx_final & genotype,2);
Tg_ang_diff(isnan(Tg_ang_diff))=[];
wrapTo360(rad2deg(circ_mean(deg2rad(Tg_ang_diff))))
circ_r(deg2rad(Tg_ang_diff))
circ_kappa(deg2rad(Tg_ang_diff))

[p,U2_obs,U2_H0]=watsons_U2_perm_test(WT_ang_diff,Tg_ang_diff,10)

% First session difference 
fig = figure;
subplot(1,2,1)
fig. Color = [1 1 1];
cdot_plot(ang_diff(HD_idx_final & ~region & ~genotype,2), 6, [.70 .70 .70],60, 10)
cdot_plot(ang_diff(HD_idx_final & region & ~genotype,2), 6, [.25 .25 .25],60, 10,gca,'-a')
subplot(1,2,2)
fig. Color = [1 1 1];
cdot_plot(ang_diff(HD_idx_final & ~region & genotype,2), 6, rgb('Orange'), 60, 10) %cortex
cdot_plot(ang_diff(HD_idx_final & region & genotype,2), 6, rgb('OrangeRed'), 60, 10,gca,'-a') %atn 

stability_idx = strcmp(groupid(:,4),'stability_rotation');

% Second vs Rotated Session 
WT_ang_diff = ang_diff(HD_idx_final & stability_idx & ~genotype,3);
WT_ang_diff(isnan(WT_ang_diff))=[];
wrapTo360(rad2deg(circ_mean(deg2rad(WT_ang_diff))))
circ_r(deg2rad(WT_ang_diff))
circ_kappa(deg2rad(WT_ang_diff))

Tg_ang_diff = ang_diff(HD_idx_final & stability_idx & genotype,3);
Tg_ang_diff(isnan(Tg_ang_diff))=[];
wrapTo360(rad2deg(circ_mean(deg2rad(Tg_ang_diff))))
circ_r(deg2rad(Tg_ang_diff))
circ_kappa(deg2rad(Tg_ang_diff))
[p,U2_obs,U2_H0]=watsons_U2_perm_test(WT_ang_diff,Tg_ang_diff,10)
fig = figure;
subplot(1,2,1)
fig. Color = [1 1 1];
cdot_plot(ang_diff(HD_idx_final & stability_idx & ~region & ~genotype,3), 6, [.70 .70 .70], 60, 10)
cdot_plot(ang_diff(HD_idx_final & stability_idx & region & ~genotype,3), 6, [.25 .25 .25], 60, 10,gca,'-a')
subplot(1,2,2)
cdot_plot(ang_diff(HD_idx_final & stability_idx & ~region & genotype,3), 6, rgb('Orange'), 60, 10)
cdot_plot(ang_diff(HD_idx_final & stability_idx & region & genotype,3), 6, rgb('OrangeRed'), 60, 10,gca,'-a')


% Rotated vs Standard 3 
WT_ang_diff = ang_diff(HD_idx_final & stability_idx & ~genotype,4);
WT_ang_diff(isnan(WT_ang_diff))=[];
wrapTo360(rad2deg(circ_mean(deg2rad(WT_ang_diff))))
circ_r(deg2rad(WT_ang_diff))
circ_kappa(deg2rad(WT_ang_diff))

Tg_ang_diff = ang_diff(HD_idx_final & stability_idx & genotype,4);
Tg_ang_diff(isnan(Tg_ang_diff))=[];
wrapTo360(rad2deg(circ_mean(deg2rad(Tg_ang_diff))))
circ_r(deg2rad(Tg_ang_diff))
circ_kappa(deg2rad(Tg_ang_diff))

[p,U2_obs,U2_H0]=watsons_U2_perm_test(WT_ang_diff,Tg_ang_diff,10)

fig = figure;
subplot(1,2,1)
fig. Color = [1 1 1];
cdot_plot(ang_diff(HD_idx_final & stability_idx & ~region & ~genotype,4), 6, [.70 .70 .70], 60, 10)
cdot_plot(ang_diff(HD_idx_final & stability_idx & region & ~genotype,4), 6, [.25 .25 .25], 60, 10,gca,'-a')
subplot(1,2,2)
cdot_plot(ang_diff(HD_idx_final & stability_idx & ~region & genotype,4), 6, rgb('Orange'), 60, 10)
cdot_plot(ang_diff(HD_idx_final & stability_idx & region & genotype,4), 6, rgb('OrangeRed'), 60, 10,gca,'-a')

% Standard 1 vs Standard 3
WT_ang_diff = ang_diff(HD_idx_final & stability_idx & ~genotype,5);
WT_ang_diff(isnan(WT_ang_diff))=[];
wrapTo360(rad2deg(circ_mean(deg2rad(WT_ang_diff))))
circ_r(deg2rad(WT_ang_diff))
circ_kappa(deg2rad(WT_ang_diff))

Tg_ang_diff = ang_diff(HD_idx_final & stability_idx & genotype,5);
Tg_ang_diff(isnan(Tg_ang_diff))=[];
wrapTo360(rad2deg(circ_mean(deg2rad(Tg_ang_diff))))
circ_r(deg2rad(Tg_ang_diff))
circ_kappa(deg2rad(Tg_ang_diff))

[p,U2_obs,U2_H0]=watsons_U2_perm_test(WT_ang_diff,Tg_ang_diff,10)

fig = figure;
subplot(1,2,1)
fig. Color = [1 1 1];
cdot_plot(ang_diff(HD_idx_final & stability_idx & ~region & ~genotype,5), 6, [.70 .70 .70], 60, 10)
cdot_plot(ang_diff(HD_idx_final & stability_idx & region & ~genotype,5), 6, [.25 .25 .25], 60, 10,gca,'-a')
subplot(1,2,2)
cdot_plot(ang_diff(HD_idx_final & stability_idx & ~region & genotype,5), 6, rgb('Orange'), 60, 10)
cdot_plot(ang_diff(HD_idx_final & stability_idx & region & genotype,5), 6, rgb('OrangeRed'), 60, 10,gca,'-a')


    specific_idx = strcmp(groupid(:,1),'LB07_S20190410173638.mat') & ...
        dir_measures_all(:,contains(vars,'mean_vector_length'),1) > .2;  
    
    stability_idx = strcmp(groupid(:,4),'dim_cue_rotation') & dir_measures_all(:,contains(vars,'mean_vector_length'),1) > .2;
    temp_id = groupid(stability_idx & genotype,:);
    
%     temp_id = groupid(specific_idx,:);
    
    for ii = 1:length(temp_id)
        temp=load(['F:\ClarkP30_Recordings\ProcessedData\',temp_id{ii,1}]);
        
        cell=find(contains(temp.spikesID.TetrodeNum,temp_id{ii,2}) & ismember(temp.spikesID.CellNum,str2double(temp_id{ii,3})))';
        
        ses=size(temp.events,2);
        theta = 0:.01:2*pi;
        color=hsv(length(theta));
    
        for iii = 1:ses
            [data_video_spk,~]=createframes_w_spikebinary(temp,iii,cell);
            fig=figure('Name',[temp.rat,'  ',temp.sessionID,'  ',...
                temp_id{ii,2},' Cell: ',num2str(temp.spikesID.CellNum(cell))],'NumberTitle','off');
            fig.Color = [1 1 1];
            subaxis(3,1,1)
            postprocessFigures.plot_HD_tuning(temp,iii,cell)
            
            subaxis(3,1,2)
            ts=data_video_spk(:,1);
            y=data_video_spk(:,3);
            x=data_video_spk(:,2);
            spkbinary=logical(data_video_spk(:,6));
            plot(x,y,'.k');
            axis image
            hold on;box off; axis off
            scatter(x(spkbinary),y(spkbinary),20,...
                interp1(rad2deg(theta)',color,data_video_spk(spkbinary,4)),'filled');
            set(gca,'FontSize',16)
            
            subaxis(3,1,3)
            max_lag = 0.05;
            % t_bin=0.005;
            t_bin=0.002;
            % Acor - taken from intrinsic frequency 2
            if t_bin / mod(max_lag, t_bin) ~= 2 % set lags so it is 'even' (odd number of coefficients and zero centered')
                max_lag = t_bin*floor(max_lag/t_bin)+max_lag*t_bin;
            end
            ses_idx = temp.Spikes{cell,1} > temp.events(1,iii) & temp.Spikes{cell,1} < temp.events(2,iii);  
            [cor, lag] = CrossCorr(temp.Spikes{cell,1}(ses_idx),temp.Spikes{cell,1}(ses_idx),'lag',...
                [-max_lag max_lag],'binsize',t_bin,'norm','prob');
            bar(lag,cor,'k')
            xlabel('Time (ms)')
            box off
            set(gcf,'OuterPosition',[1856 361 338 605])
            %             imAlpha=ones(size(temp.ratemap{cell,iii}));
            %             imAlpha(isnan(temp.ratemap{cell,iii}))=0;
            %             imagesc(temp.ratemap{cell,iii},'AlphaData',imAlpha);
            %             axis xy; axis off; hold on; box off; axis image;
            %             colormap(gca,viridis(255))
            %
            %
        end
    end
    
    



% BD_cell_data = dir_measures_all(HD_idx & ~place_idx & BD_idx,:,:);
% BD_cell_id = groupid(HD_idx & ~place_idx & BD_idx,:);
% angBins=0:6:360;
% bin_centers=movmedian(angBins,2);
% bin_centers(1)=[];
% for i = 1:length(HD_cell_id)
% 
%     temp=load(['F:\ClarkP30_Recordings\ProcessedData\',HD_cell_id{i,1}]);
% 
%     cell=find(contains(temp.spikesID.TetrodeNum,HD_cell_id{i,2}) & ismember(temp.spikesID.CellNum,str2double(HD_cell_id{i,3})))';
% 
%     ses=size(temp.events,2);
% 
%     [data_video_spk,~]=createframes_w_spikebinary(temp,1,cell);
%     [~,within,~]=HD_cell_analysis.four_quarter_stability(data_video_spk,temp.samplerate,'std');
% 
%     max_lag = 0.05;
%     % t_bin=0.005;
%     t_bin=0.002;
%     % Acor - taken from intrinsic frequency 2
%     if t_bin / mod(max_lag, t_bin) ~= 2 % set lags so it is 'even' (odd number of coefficients and zero centered')
%         max_lag = t_bin*floor(max_lag/t_bin)+max_lag*t_bin;
%     end
%     [cor, lag] = CrossCorr(temp.Spikes{cell,1},temp.Spikes{cell,1},'lag',[-max_lag max_lag],'binsize',t_bin,'norm','prob');
% 
%       % Plot Spikes on Path by HD
%     fig = figure;
%     fig.Color = [1 1 1];
%     theta = 0:.01:2*pi;
%     color=hsv(length(theta));
% 
%     subaxis(1,2,1)
%     ts=data_video_spk(:,1);
%     y=data_video_spk(:,3);
%     x=data_video_spk(:,2);
%     spkbinary=logical(data_video_spk(:,6));
%     plot(x,y,'.k');
%     axis image
%     hold on;box off; axis off
%     scatter(x(spkbinary),y(spkbinary),20,...
%         interp1(rad2deg(theta)',color,data_video_spk(spkbinary,4)),'filled');
%     title([extractBefore(HD_cell_id{i,1},'_'),' cell: ',HD_cell_id{i,3}])
%     set(gca,'FontSize',16)
% 
%     subaxis(1,2,2)
%     postprocessFigures.plot_HD_tuning(temp,1,cell)
%     set(gca,'FontSize',16)
% 
%     fig = figure;
%     fig.Color = [1 1 1];
% 
%     subaxis(2,1,1)
%     sam=size(temp.avgwave{cell,1},2);
%     waves=zeros(4,sam);
%     waves(1:size(temp.avgwave{cell,1},1),:)=temp.avgwave{cell,1};
%     plot(1:sam,waves(1,:),'LineWidth',2, 'color','k');hold on
%     plot(sam+1:sam*2,waves(2,:),'LineWidth',2, 'color','k');
%     plot(1:sam,waves(3,:)+max(abs([waves(1,:),waves(2,:)])),'LineWidth',2, 'color','k');
%     plot(sam+1:sam*2,waves(4,:)+max(abs([waves(1,:),waves(2,:)])),'LineWidth',2, 'color','k');
%     axis tight
%     hold on;box off; axis off
%     subaxis
% 
% 
%     subaxis(2,1,2)
%     bar(lag,cor,'k')
%     xlabel('Time (ms)')
%     box off
% 
% 
% %
% %     print(fig,'-dpng', '-r600',['d:\Users\BClarkLab\Desktop\Laura Temp\WT_HD\',...
% %         HD_cell_id{i,1:3},'.png'])
% %     close
% %
% %
% end




% HD_cell_data = dir_measures_all(HD_idx & ~place_idx & ~BD_idx & genotype,:,:);
% HD_cell_id = groupid(HD_idx & ~place_idx & ~BD_idx & genotype,:);
% angBins=0:6:360;
% bin_centers=movmedian(angBins,2);
% bin_centers(1)=[];
% for i = 1:length(HD_cell_id)
%
%     temp=load(['F:\ClarkP30_Recordings\ProcessedData\',HD_cell_id{i,1}]);
%
%     cell=find(contains(temp.spikesID.TetrodeNum,HD_cell_id{i,2}) & ismember(temp.spikesID.CellNum,str2double(HD_cell_id{i,3})))';
%
%     ses=size(temp.events,2);
%
%     [data_video_spk,~]=createframes_w_spikebinary(temp,1,cell);
%     [~,within,~]=HD_cell_analysis.four_quarter_stability(data_video_spk,temp.samplerate,'std');
%
%     max_lag = 0.05;
%     % t_bin=0.005;
%     t_bin=0.002;
%     % Acor - taken from intrinsic frequency 2
%     if t_bin / mod(max_lag, t_bin) ~= 2 % set lags so it is 'even' (odd number of coefficients and zero centered')
%         max_lag = t_bin*floor(max_lag/t_bin)+max_lag*t_bin;
%     end
%     [cor, lag] = CrossCorr(temp.Spikes{cell,1},temp.Spikes{cell,1},'lag',[-max_lag max_lag],'binsize',t_bin,'norm','prob');
%
%     if HD_cell_data(i,contains(vars,'mean_vector_length'),1) >= .45
%       % Plot Spikes on Path by HD
%     fig = figure;
%     fig.Color = [1 1 1];
%     theta = 0:.01:2*pi;
%     color=hsv(length(theta));
%
%     subaxis(1,2,1)
%     ts=data_video_spk(:,1);
%     y=data_video_spk(:,3);
%     x=data_video_spk(:,2);
%     spkbinary=logical(data_video_spk(:,6));
%     plot(x,y,'.k');
%     axis image
%     hold on;box off; axis off
%     scatter(x(spkbinary),y(spkbinary),20,...
%         interp1(rad2deg(theta)',color,data_video_spk(spkbinary,4)),'filled');
%     title([extractBefore(HD_cell_id{i,1},'_'),' cell: ',HD_cell_id{i,3}])
%     set(gca,'FontSize',16)
%
%     subaxis(1,2,2)
%     postprocessFigures.plot_HD_tuning(temp,1,cell)
%     set(gca,'FontSize',16)
%
%     fig = figure;
%     fig.Color = [1 1 1];
%
%     subaxis(2,1,1)
%     sam=size(temp.avgwave{cell,1},2);
%     waves=zeros(4,sam);
%     waves(1:size(temp.avgwave{cell,1},1),:)=temp.avgwave{cell,1};
%     plot(1:sam,waves(1,:),'LineWidth',2, 'color','k');hold on
%     plot(sam+1:sam*2,waves(2,:),'LineWidth',2, 'color','k');
%     plot(1:sam,waves(3,:)+max(abs([waves(1,:),waves(2,:)])),'LineWidth',2, 'color','k');
%     plot(sam+1:sam*2,waves(4,:)+max(abs([waves(1,:),waves(2,:)])),'LineWidth',2, 'color','k');
%     axis tight
%     hold on;box off; axis off
%
%
%     subaxis(2,1,2)
%     bar(lag,cor,'k')
%     xlabel('Time (ms)')
%     box off
%
%     else
%         continue
%     end
%
% %
% %     print(fig,'-dpng', '-r600',['d:\Users\BClarkLab\Desktop\Laura Temp\WT_HD\',...
% %         HD_cell_id{i,1:3},'.png'])
% %     close
% %
% %
% end



% Place_cell_data = dir_measures_all(~HD_idx & place_idx & ~BD_idx & ~IN_idx & ~genotype,:,:);
% Place_cell_id = groupid(~HD_idx & place_idx & ~BD_idx & ~IN_idx & ~genotype,:);
% 
% for i=1:size(Place_cell_id,1)
%     
%     if  Place_cell_data(i,contains(vars,'InformationContent'),1) >= .4
%         
%         tt = str2double(regexp(Place_cell_id{i,2}, '\d*','Match'));
%         cylinder_example(Place_cell_id{i,1},tt,str2double(Place_cell_id{i,3}),1)
%         set(findall(gcf,'-property','FontSize'),'FontSize',27)
%         
%     else
%         continue
%     end
% end
% 
% 
% function cylinder_example(session,tt,cell,ns)
% data = load(session);
% [i]=find_cells(data,tt,cell);
% tetrode=strsplit(data.spikesID.paths{i},filesep);
% tetrode=tetrode{end};
% trodeID=str2double(extractBetween(tetrode,'TT','.'));
% colorcode = 'r';
% fig=figure('Name',[data.rat,'  ',data.sessionID,'  ',...
%     tetrode,' Cell: ',num2str(data.spikesID.CellNum(i))],'NumberTitle','off');
% fig.Color = [1 1 1];
% p = panel(fig);
% p.pack(4, 1);
% % set margins
% p.de.margin = 6;
% p(1, 1).select();
% ax = gca;
% [data_video_spk,data_video_nospk]=createframes_w_spikebinary(data,ns,i);
% postprocessFigures.spikesonpath_2d(ax,data_video_spk,data.lfp.ts,data.lfp.theta_phase(trodeID,:),colorcode)
% title('')
% disp(sprintf('IC: %4.2f  %4.2f hz',...
%     data.measures(contains(data.varnames,["InformationContent","PeakRate"]))))
% p(2, 1).select();
% ax = gca;
% postprocessFigures.ratemaps_2d(ax,data.ratemap{i,ns},data.measures(i,:,ns),data.varnames)
% title('')
% binside=mean([range(data_video_nospk(:,2))/length(data.ratemap{i,ns}),...
%     range(data_video_nospk(:,3))/length(data.ratemap{i,ns})]);
% results=pass_index(data_video_nospk(:,1),data_video_nospk(:,2:3),...
%     data_video_spk(data_video_spk(:,6)==1,1),...
%     [data.lfp.ts(data.lfp.ts>=data.events(1,ns) & data.lfp.ts<=data.events(2,ns))]',...
%     [data.lfp.signal(trodeID,data.lfp.ts>=data.events(1,ns) & data.lfp.ts<=data.events(2,ns))]',...
%     'plots',0,'method','place','binside',round(binside),'sample_along','arc_length');
% bins=length(data.ratemap{i,ns});
% xedge=linspace(-1,1,bins+1);
% phaseedge=linspace(0,720,bins*4);
% phase_spk = interp1(data.lfp.ts,data.lfp.theta_phase(trodeID,:),data_video_spk(data_video_spk(:,6)==1,1));
% x_spk = interp1(results.ts,results.pass_index,data_video_spk(data_video_spk(:,6)==1,1));
% spkmap=histcounts2([x_spk;x_spk],[phase_spk;phase_spk+2*pi]*180/pi,xedge,phaseedge);
% phase_occ = interp1(data.lfp.ts,data.lfp.theta_phase(trodeID,:),data_video_spk(data_video_spk(:,6)==0,1));
% x_occ = interp1(results.ts,results.pass_index,data_video_spk(data_video_spk(:,6)==0,1));
% occ=histcounts2([x_occ;x_occ],[phase_occ;phase_occ+2*pi]*180/pi,xedge,phaseedge);
% occ = occ/data.samplerate;
% phasemap=spkmap./occ;
% phasemap(isnan(phasemap)) = 0;
% phasemap(isinf(phasemap)) = 0;
% h=4;
% phase_size = size(phasemap,2);
% myfilter = fspecial('gaussian',[4 24]*h, h);
% phasemap = imfilter([phasemap,phasemap,phasemap],myfilter,'replicate');
% phasemap = phasemap(:,phase_size:phase_size*2);
% p(3, 1).select();
% scatter([x_spk;x_spk],[phase_spk;phase_spk+2*pi],20,'Filled','k');
% axis tight
% axis square
% axis off
% p(4, 1).select();
% pcolor(flipud(rot90(phasemap)));shading flat;box off;axis off;axis tight;
% colormap(viridis(255))
% axis square
% axis on
% ylabel('\theta Phase')
% xlabel('')
% % set(gca,'XTick',linspace(ax.XLim(1),ax.XLim(2),3),'XTickLabels',[-1,0,1])
% set(gca,'XTickLabels',[])
% ax = gca;
% set(gca,'YTick',linspace(ax.YLim(1),ax.YLim(2),3),'YTickLabels',{'0','2\pi','4\pi'})
% % set(gcf,'OuterPosition',[1050 34 369 1017])?
% end
% 
% 
% 
% fig = figure; 
% fig.Color = [1 1 1]
% stability_idx = strcmp(groupid(:,4),'stability') | strcmp(groupid(:,4),'stability_rotation');
% scatter(dir_measures_all(HD_idx_final & stability_idx & ~genotype,contains(vars,'preferred_Direction'),1),...
%     dir_measures_all(HD_idx_final & stability_idx & ~genotype,contains(vars,'preferred_Direction'),2),100,'filled','k')
% hold on
% scatter(dir_measures_all(HD_idx_final & stability_idx & genotype,contains(vars,'preferred_Direction'),1),...
%     dir_measures_all(HD_idx_final & stability_idx & genotype,contains(vars,'preferred_Direction'),2),100,'filled','r')
% hold on; 
% plot(1:1:360,1:1:360,'--k','LineWidth',2)
% xlabel('Standard 1')
% ylabel('Standard 2')
% title('Preferred Direction between Standard Sessions')
% legend({'F344','TgF344-AD'})
% legend boxoff
% set(gca,'LineWidth',2,'FontSize',16,'FontName','Helvetica','FontWeight','Bold')
% 
% fig = figure; 
% fig.Color = [1 1 1]
% stability_idx = strcmp(groupid(:,4),'stability_rotation');
% scatter(dir_measures_all(HD_idx_final & stability_idx & ~genotype,contains(vars,'preferred_Direction'),2),...
%     dir_measures_all(HD_idx_final & stability_idx & ~genotype,contains(vars,'preferred_Direction'),3),200,'filled',...
%     'MarkerEdgeColor',[.25 .25 .25],'MarkerFaceColor',[.25 .25 .25])
% hold on
% scatter(dir_measures_all(HD_idx_final & stability_idx & genotype,contains(vars,'preferred_Direction'),2),...
%     dir_measures_all(HD_idx_final & stability_idx & genotype,contains(vars,'preferred_Direction'),3),200,'filled',...
%     'MarkerEdgeColor',rgb('DarkOrange'),'MarkerFaceColor',rgb('DarkOrange'))
% hold on; 
% plot([0 270],[90 360],'--k','LineWidth',2)
% xlim([0 361])
% ylim([0 361])
% plot([90 360],[0 270],'--k','LineWidth',2)
% legend({' ',' ',' '})
% legend boxoff
% set(gca,'LineWidth',2,'FontSize',16,'FontName','Helvetica','FontWeight','Bold','OuterPosition',[0 0 1 1])
% 
% 
% 
% fig = figure; 
% fig.Color = [1 1 1]
% stability_idx = strcmp(groupid(:,4),'stability_rotation');
% scatter(dir_measures_all(HD_idx_final & stability_idx & ~genotype,contains(vars,'preferred_Direction'),1),...
%     dir_measures_all(HD_idx_final & stability_idx & ~genotype,contains(vars,'preferred_Direction'),2),200,'filled',...
%     'MarkerEdgeColor',[.25 .25 .25],'MarkerFaceColor',[.25 .25 .25])
% hold on
% scatter(dir_measures_all(HD_idx_final & stability_idx & genotype,contains(vars,'preferred_Direction'),1),...
%     dir_measures_all(HD_idx_final & stability_idx & genotype,contains(vars,'preferred_Direction'),2),200,'filled',...
%     'MarkerEdgeColor',rgb('DarkOrange'),'MarkerFaceColor',rgb('DarkOrange'))
% hold on; 
% plot(1:1:360,1:1:360,'--k','LineWidth',2)
% xlim([0 361])
% ylim([0 361])
% legend({' ',' ',' '})
% legend boxoff
% set(gca,'LineWidth',2,'FontSize',16,'FontName','Helvetica','FontWeight','Bold','OuterPosition',[0 0 1 1])
% 



