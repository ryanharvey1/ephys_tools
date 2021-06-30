% manual_ripple_classification
basepath = 'F:\Projects\PAE_PlaceCell';

if exist(fullfile(basepath,'SWR_workspace.mat'),'file')
    load(fullfile(basepath,'SWR_workspace.mat'),'SWR');
else
    sessions = dir(fullfile(basepath,'SWR','*.mat'));
    maps = [];
    maps_filtered = [];
    duration = [];
    peakFrequency = [];
    peakAmplitude = [];
    session_id = [];
    ripple_id = [];
    for s = 1:length(sessions)
        disp(fullfile(sessions(s).folder,sessions(s).name))
        ripple_info = load(fullfile(sessions(s).folder,sessions(s).name));
        
        % check if ripple maps are 151 in length (there's a bug somewhere in Sync)
        if size(ripple_info.maps.ripples,2) ~= 151
            ripple_info = sync_to_size(ripple_info);
        end
        
        maps = [maps;ripple_info.maps.unfiltered_ripples];
        maps_filtered = [maps_filtered;ripple_info.maps.ripples];
        
        duration = [duration;ripple_info.ripple_data.duration];
        peakFrequency = [peakFrequency;ripple_info.ripple_data.peakFrequency];
        peakAmplitude = [peakAmplitude;ripple_info.ripple_data.peakAmplitude];
        
        session_id = [session_id; repmat({ripple_info.ripples.detectorinfo.ProcessedDatafile},...
            size(ripple_info.maps.ripples,1),1)];
        
        ripple_id = [ripple_id;[1:size(ripple_info.maps.ripples,1)]'];
    end
    SWR.maps = maps;
    SWR.maps_zscore = zscore(SWR.maps,[],2);
    SWR.maps_filtered = maps_filtered;
    
    SWR.data.duration = duration;
    SWR.data.peakFrequency = peakFrequency;
    SWR.data.peakAmplitude = peakAmplitude;
    SWR.data.session_id = session_id;
    SWR.data.ripple_id = ripple_id;
    SWR.data.result = zeros(size(SWR.maps,1),1);
    SWR.data.labeled = zeros(size(SWR.maps,1),1);
    
    SWR.data = struct2table(SWR.data);
    
    save(fullfile(basepath,'SWR_workspace.mat'),'SWR');
end
%%
move_fig(basepath)

%%
diff(SWR.maps_zscore,1,2);

 [coeff,score,latent,tsquared,explained,mu] = pca(normalize(SWR.maps,2,'center'));
   figure
    h = scatter3(score(:,1),score(:,2),score(:,3));

%%
window_ = round(size(SWR.maps,2)/2-25:size(SWR.maps,2)/2+25);

minPCvar = 0.8;
probClus = 0;
autoClus = 1;
numRep = 100;
waveletParams = {1000, [70 300], 1};
save_dir = 'F:\Projects\PAE_PlaceCell\som_ripple_figs';
[clusData, clussMap, clusNum, probMatrix] = RhythSOM_Classifier(normalize(SWR.maps(:,window_),2,'range'), minPCvar,...
    probClus, autoClus, numRep, waveletParams,save_dir) ;

%%
list_of_good = find(clusData == 28);
figure
imagesc(SWR.maps_zscore(list_of_good,:))
figure;
plot(SWR.maps_zscore(list_of_good,:)')


list_of_good = find(clusData == 4);
for i = 1:length(list_of_good)
    fig =  figure;
    plot(SWR.maps(list_of_good(i),:))
    pause
    close(fig)
end

%%

% figure;
% 
% [coeff,score,latent,tsquared,explained,mu] = pca(SWR.maps_zscore);
% subplot(3,2,5)
% h = scatter(score(:,1),score(:,2),3,zeros(length(score(:,2)),3)+.5,'Filled');
% h.SizeData = zeros(length(score(:,2)),1) + 5;
% h.MarkerFaceAlpha = .7;
% axis tight
% hold on
% 
% for i = 1:length(SWR.maps)
%     
%     h.CData = zeros(length(score(:,2)),3);
%     
%     sgtitle(['duration: ',num2str(SWR.duration(i)),'s ',...
%         'freq: ',num2str(round(SWR.peakFrequency(i))),'Hz ',...
%         num2str(i),' of ',num2str(size(SWR.maps,1))],...
%         'color','r')
%     
%     subplot(3,2,5)
%     h2 = scatter(score(i,1),score(i,2),10,[1 0 0],'Filled');
%     
%     
%     subplot(3,2,1)
%     plot(1:151,SWR.maps_zscore(i,:),'Color',[.7,.7,.7])
%     
%     subplot(3,2,2)
%     plot(1:151,SWR.maps_filtered(i,:),'Color',[.7,.7,.7])
%     darkBackground(gcf,[0.1 0.1 0.1],[0.7 0.7 0.7])
%     
%     SWR.data.result(i) = input('ripple ? ');
%     SWR.data.labeled(i) = 1;
%     delete(h2)
%     
%     if SWR.result(i)
%         subplot(3,2,3)
%         plot(1:151,SWR.maps_zscore(i,:),...
%             'Color',[rand(1, 1),rand(1, 1),rand(1, 1),0.1])
%         title('good')
%         hold on
%     end
%     if ~SWR.result(i)
%         subplot(3,2,4)
%         plot(1:151,SWR.maps_zscore(i,:),...
%             'Color',[rand(1, 1),rand(1, 1),rand(1, 1),0.1])
%         title('bad')
%         hold on
%     end
%     
%     
%     darkBackground(gcf,[0.1 0.1 0.1],[0.7 0.7 0.7])
% end



% maps_zscore = zscore(maps,[],2);
% 
% 
% [reduction, umap, clusterIdentifiers] = run_umap(maps_zscore);
% clusters = unique(clusterIdentifiers);
% p = ceil(sqrt(length(clusters)));
% figure;
% for i = double(clusters)
%     subplot(p,p,i+1)
%     plot(1:151,maps(clusterIdentifiers == i,:),'Color',[.7,.7,.7,0.1])
%     box off
%     axis off
% end
% darkBackground(gcf,[0.1 0.1 0.1],[0.7 0.7 0.7])
% 
% figure;
% for i = double(clusters)
%     plot(1:151,maps_zscore(clusterIdentifiers == i,:),...
%         'Color',[rand(1, 1),rand(1, 1),rand(1, 1),0.1])
%     hold on
% end
% darkBackground(gcf,[0.1 0.1 0.1],[0.7 0.7 0.7])
% axis off
% 
% 
% idx = kmeans(reduction,3,'Distance','cityblock');
% figure;
% scatter(reduction(idx==1,1),reduction(idx==1,2),'r','Filled')
% hold on
% scatter(reduction(idx==2,1),reduction(idx==2,2),'g','Filled')
% scatter(reduction(idx==3,1),reduction(idx==3,2),'b','Filled')
% axis image
% 
% 
% 
% 
% figure;
% scatter3(duration,peakFrequency,peakAmplitude,'k','Filled')
% 
% maps_filtered_zscore = zscore(maps_filtered,[],2);
% 
% 
% maps_zscore = zscore(maps,[],2);
% 
% [coeff,score,latent,tsquared,explained,mu] = pca(maps_zscore);
% [coeff,score,latent,tsquared,explained,mu] = pca(maps);
% 
% 
% figure;
% scatter3(score(:,1),score(:,2),score(:,3),'k','Filled')
% axis image
% 
% idx = kmeans(score,2);
% figure;
% scatter3(score(:,1),score(:,2),score(:,3),'k','Filled')
% hold on
% scatter3(score(idx==1,1),score(idx==1,2),score(idx==1,3),'r','Filled')
% axis image
% 
% figure
% 
% biplot(coeff(:,1:3),'scores',score(:,1:3));
% axis image

function ripple_info = sync_to_size(ripple_info)
for i = 1:size(ripple_info.maps.ripples,1)
    
    ripples(i,:) = interp1(linspace(1,151,size(ripple_info.maps.ripples,2)),...
        ripple_info.maps.ripples(i,:),1:151);
    
    frequency(i,:) = interp1(linspace(1,151,size(ripple_info.maps.frequency,2)),...
        ripple_info.maps.frequency(i,:),1:151);
    
    phase(i,:) = interp1(linspace(1,151,size(ripple_info.maps.phase,2)),...
        ripple_info.maps.phase(i,:),1:151);
    
    amplitude(i,:) = interp1(linspace(1,151,size(ripple_info.maps.amplitude,2)),...
        ripple_info.maps.amplitude(i,:),1:151);
    
    unfiltered_ripples(i,:) = interp1(linspace(1,151,size(ripple_info.maps.unfiltered_ripples,2)),...
        ripple_info.maps.unfiltered_ripples(i,:),1:151);
end
ripple_info.maps.ripples = ripples;
ripple_info.maps.frequency = frequency;
ripple_info.maps.phase = phase;
ripple_info.maps.amplitude = amplitude;
ripple_info.maps.unfiltered_ripples = unfiltered_ripples;
end

function [] = move_fig(basepath)
% move figure with arrow keys.
S.fh = figure('units','pixels',...
    'menubar','none',...
    'name','move_fig',...
    'numbertitle','off',...
    'keypressfcn',@fh_kpfcn);
S.basepath = basepath;
S.start_of_session = 1;
S.i = 1;
% S.tx = uicontrol('style','text',...
%     'units','pixels',...
%     'position',[60 120 80 20],...
%     'fontweight','bold');
guidata(S.fh,S)
disp('Hit Enter')
end
function [S] = fh_kpfcn(H,E)
% Figure keypressfcn
S = guidata(H);
P = get(S.fh,'position');

i = S.i;

if S.start_of_session
    load(fullfile(S.basepath,'SWR_workspace.mat'),'SWR');
    
    i = find(SWR.data.labeled,1,'last');
    if isempty(i)
        i = 1;
    end
    S.SWR = SWR;
else
    SWR = S.SWR;
end

if ~S.start_of_session
    switch E.Key
        case 'rightarrow'
            i = i + 1;
        case 'leftarrow'
            i = i - 1;
        case 'uparrow'
            SWR.data.result(i) = 1;
            SWR.data.labeled(i) = 1;
            i = i + 1;
        case 'downarrow'
            SWR.data.result(i) = 0;
            SWR.data.labeled(i) = 1;
            i = i + 1;
        otherwise
            save(fullfile(S.basepath,'SWR_workspace.mat'),'SWR');
    end
end
disp([SWR.data.session_id{i},...
    ' ripple ',num2str(SWR.data.ripple_id(i)),...
    ' Labeled: ',num2str(SWR.data.labeled(i)),...
    ' ripple: ',num2str(SWR.data.result(i))])

if S.start_of_session
    [coeff,score,latent,tsquared,explained,mu] = pca(SWR.maps);
    subplot(3,2,5)
    h = scatter(score(:,1),score(:,2),3,zeros(length(score(:,2)),3)+.5,'Filled');
    h.SizeData = zeros(length(score(:,2)),1) + 5;
    h.MarkerFaceAlpha = .7;
    axis tight
    hold on
    h.CData = zeros(length(score(:,2)),3);
    S.score = score;
    xlabel('pc1')
    ylabel('pc2')
    
    subplot(3,2,6)
    h = scatter(SWR.data.duration,SWR.data.peakFrequency,3,...
        zeros(length(score(:,2)),3)+.5,'Filled');
    h.SizeData = zeros(length(score(:,2)),1) + 5;
    h.MarkerFaceAlpha = .7;
    axis tight
    hold on
    h.CData = zeros(length(score(:,2)),3);
    xlabel('duration(s)')
    ylabel('freq')
end
S.start_of_session = 0;

sgtitle(['duration: ',num2str(SWR.data.duration(i)),'s ',...
    'freq: ',num2str(round(SWR.data.peakFrequency(i))),'Hz ',...
    num2str(i),' of ',num2str(size(SWR.maps,1))],...
    'color','r')

subplot(3,2,1)
plot(1:151,SWR.maps_zscore(i,:),'Color',[.7,.7,.7])
title('raw ripple')

subplot(3,2,2)
plot(1:151,SWR.maps_filtered(i,:),'Color',[.7,.7,.7])
title('filtered ripple')

subplot(3,2,5)
scatter(S.score(i,1),S.score(i,2),10,[1,1,1],'Filled');

subplot(3,2,6)
scatter(SWR.data.duration(i),SWR.data.peakFrequency(i),10,[1,1,1],'Filled');

darkBackground(gcf,[0.1 0.1 0.1],[0.7 0.7 0.7])

if i > 1
    if SWR.data.result(i-1) && SWR.data.labeled(i-1)
        subplot(3,2,3)
        plot(1:151,SWR.maps_zscore(i-1,:),...
            'Color',[rand(1, 1),rand(1, 1),rand(1, 1),0.1])
        title('good','Color','g')
        hold on
        
        subplot(3,2,5)
        scatter(S.score(i-1,1),S.score(i-1,2),10,[0 1 0],'Filled');
    end
    if ~SWR.data.result(i-1) && SWR.data.labeled(i-1)
        subplot(3,2,4)
        plot(1:151,SWR.maps_zscore(i-1,:),...
            'Color',[rand(1, 1),rand(1, 1),rand(1, 1),0.1])
        title('bad','Color','r')
        hold on
        
        subplot(3,2,5)
        scatter(S.score(i-1,1),S.score(i-1,2),10,[1 0 0],'Filled');
    end
 
    if SWR.data.result(i-1) && SWR.data.labeled(i-1)
        subplot(3,2,6)
        scatter(SWR.data.duration(i-1),SWR.data.peakFrequency(i-1),10,[0 1 0],'Filled');
    end
    
    if ~SWR.data.result(i-1) && SWR.data.labeled(i-1)
        subplot(3,2,6)
        scatter(SWR.data.duration(i-1),SWR.data.peakFrequency(i-1),10,[1 0 0],'Filled');
    end
end
darkBackground(gcf,[0.1 0.1 0.1],[0.7 0.7 0.7])

S.i = i;
S.SWR = SWR;
guidata(S.fh,S)

end
