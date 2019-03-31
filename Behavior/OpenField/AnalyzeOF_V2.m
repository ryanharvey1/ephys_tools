% contains(param_idx,'lgOF') & contains(param_idx,'Tg');
% 
% dwellQuad_tg1=vertcat(params.dwellQuad{contains(param_idx,'lgOF') & contains(param_idx,'Tg') & contains(param_idx,'day1')});
% dwellQuad_wt1=vertcat(params.dwellQuad{contains(param_idx,'lgOF') & contains(param_idx,'Wt') & contains(param_idx,'day1')});
% 
% dwellQuad_tg2=vertcat(params.dwellQuad{contains(param_idx,'lgOF') & contains(param_idx,'Tg') & contains(param_idx,'day2')});
% dwellQuad_wt2=vertcat(params.dwellQuad{contains(param_idx,'lgOF') & contains(param_idx,'Wt') & contains(param_idx,'day2')});
% 
% %Histograms for quadrant dwell time
% bincenter=[0:360/16:360]-22.5;
% bincenter=bincenter(2:end);
% %%
% f=figure; f.Color=[1 1 1];
% subplot(2,2,1)
% bar(sum(dwellQuad_tg1,1)/sum(dwellQuad_tg1(:)),'FaceColor','r');
% title('tg day 1')
% ylim([0 .6])
% xlabel('Quadrant Angle')
% ylabel('Normalized Dwell Time')
% set(gca,'XTick',[1:16],'XTickLabel',bincenter)
% 
% box off 
% 
% subplot(2,2,3)
% bar(sum(dwellQuad_wt1,1)/sum(dwellQuad_wt1(:)),'FaceColor',[.5 .5 .5]);
% title('wt day 1')
% ylim([0 .6])
% xlabel('Quadrant Angle')
% ylabel('Normalized Dwell Time')
% set(gca,'XTick',[1:16],'XTickLabel',bincenter)
% 
% box off 
% 
% subplot(2,2,2)
% bar(sum(dwellQuad_tg2,1)/sum(dwellQuad_tg2(:)),'FaceColor','r');
% title('tg day 2')
% ylim([0 .6])
% xlabel('Quadrant Angle')
% ylabel('Normalized Dwell Time')
% set(gca,'XTick',[1:16],'XTickLabel',bincenter)
% 
% box off 
% 
% subplot(2,2,4)
% bar(sum(dwellQuad_wt2,1)/sum(dwellQuad_wt2(:)),'FaceColor',[.5 .5 .5]);
% title('wt day 2')
% ylim([0 .6])
% xlabel('Quadrant Angle')
% ylabel('Normalized Dwell Time')
% set(gca,'XTick',[1:16],'XTickLabel',bincenter)
% 
% box off 
% %%
% 
% 
% figure; plot(cos(0:2*pi/1000:2*pi)*101,sin(0:2*pi/1000:2*pi)*101,'-k');
% 
% for i=1:length(y)
% hold on; plot(x{:,i},y{:,i},'--r')
% end
% 
% hold on; plot(params.transcoords{j}(:,1),params.transcoords{j}(:,2),'.b')
% 
% hold on; comet(params.transcoords{j}(:,1),params.transcoords{j}(:,2))
% 
% 
% 
% 
% 
% %%%Create a figure for Stop v Motion
% 
% 
%     figure;
%      
%        plot(params.transcoords{j}(:,1),params.transcoords{j}(:,2),'.k');hold on
%        for ii=1:length(start)
%            motionless{ii}=[params.transcoords{j}(start(ii):ends(ii),1),params.transcoords{j}(start(ii):ends(ii),2)];
%            plot(params.transcoords{j}(start(ii):ends(ii),1),params.transcoords{j}(start(ii):ends(ii),2),'.r');hold on
%        end
%            
%             title('Path with Stops')
%             ylim([0 .6])
%             xlabel('Quadrant Angle')
%             ylabel('Normalized Dwell Time')
%             set(gca,'XTick',[1:16],'XTickLabel',bincenter)
% 
% box off 
%      
        
%%"BEE SWARM" PLOT

%%Path Length 
param_idx=params.PathName;
contains(param_idx,'lgOF') & contains(param_idx,'Tg');
addpath('D:\Users\BClarkLab\GoogleDrive\MatlabDir\BClarkToolbox\Analysis\Visualize\plotSpread')


pathL_tg1=vertcat(params.pathL{contains(param_idx,'lgOF') & contains(param_idx,'Tg') & contains(param_idx,'day1')});
pathL_wt1=vertcat(params.pathL{contains(param_idx,'lgOF') & contains(param_idx,'Wt') & contains(param_idx,'day1')});

pathL_tg2=vertcat(params.pathL{contains(param_idx,'lgOF') & contains(param_idx,'Tg') & contains(param_idx,'day2')});
pathL_wt2=vertcat(params.pathL{contains(param_idx,'lgOF') & contains(param_idx,'Wt') & contains(param_idx,'day2')});

%%Varargins for the function and to customize the graph
%%Important Note: shownMW creates a line for either mean median or both, '1' creates for
%%both(used in line 116) 
data=[pathL_tg1; pathL_wt1;pathL_tg2;pathL_wt2];
group=[ones(1,12) (ones(1,12)+1)]';
catIdx=[group; group];
dayIdx = ['Day1';'Day2']; %%labels x axis
xValues_idx=[1;2];%the values of the x axis
distBeesIdx=[zeros(24,1); ones(24,1)];
yLabel_idx= ['Path Length Score']; 

figure; h=plotSpread(data,'categoryIdx',catIdx,'xNames',dayIdx,'xValues',xValues_idx,...
    'categoryColors',{'r',[.5,.5,.5]},'DistributionIdx',distBeesIdx,'yLabel',yLabel_idx,'showMM',1,...
    'distributionMarkers','o');
title('Path Length')
legend('tg','wt','mean','median')

box off



%%plot Number of Stops
for j=1:length(param_idx)
params.NumStops{j}= size(params.stops{j});
params.NumStops{j}=[params.NumStops{j}(:,2)];
end

%% Ratemaps 

for i=1:length(params.rateMap)
    if  isempty(params.rateMap{i})
        continue
    end
    
    map=params.rateMap{i};
    
    [~,fig]=PerfectCircRateMap(map,1);
    title(param_idx(i))
    
    print(gcf,'-dpng', '-r300',['d:\Users\BClarkLab\Desktop\RateMap\',erase(param_idx{i},["D:\Maria\OF_Excel\",".xlsx"]),'.png'])
    close all
    
end

%% Plot Paths

for i=1:length(params.rateMap)
    
    if  isempty(params.rateMap{i})
        continue
    end
    
    figure; plot(cos(0:2*pi/1000:2*pi)*101,sin(0:2*pi/1000:2*pi)*101,'-k');
    
    hold on; plot(params.transcoords{j}(:,1),params.transcoords{j}(:,2),'.k')
    
  %loop through home base and plot
    
%     print(gcf,'-dpng', '-r300',['d:\Users\BClarkLab\Desktop\RateMap\',erase(param_idx{i},["D:\Maria\OF_Excel\",".xlsx"]),'.png'])
%     close all
    
end