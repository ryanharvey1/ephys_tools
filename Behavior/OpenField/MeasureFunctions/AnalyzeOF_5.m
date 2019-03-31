
%%
%%%%%%%%%%% Plot Paths with home base outline and cue location
load('params_V8');
param_idx=params.PathName;
params.cueCenter{1}=[];

c=[params.transcue{19}(4,:);params.transcue{19}(2,:);params.transcue{19}(1,:);params.transcue{19}(3,:)];
cue = polyshape(c);
[cueX,cueY]=centroid(cue);
     

for j=1:length(params.transcoords)
    
    if  isempty(params.rateMap{j})
        continue
    end
    
    figure; plot(cos(0:2*pi/1000:2*pi)*101,sin(0:2*pi/1000:2*pi)*101,'-k','LineWidth',2);
    
    hold on; plot(params.transcoords{j}(:,1),params.transcoords{j}(:,2),'--k','LineWidth',2)
    
    for hb=1:size(params.HBcoords{j},1)
        hold on; plot(params.HBcoords{j}{hb,1},params.HBcoords{j}{hb,2},'r','LineWidth',2)
    end
    
    
    if contains(param_idx{j},'day2_lgOF')
        hold on; plot(cueX,cueY,'*r', 'MarkerSize',25,'LineWidth',3) 
    else
        hold on; o=plot(cue);
        o.FaceColor=[0 0 0];
        o.FaceAlpha=1;
    end
    
    axis image
    axis off
    title(erase(param_idx{j},["D:\Maria\OF_Excel\",".xlsx","_"]))
    
    print(gcf,'-dpng', '-r300',['D:\Maria\Figures\',erase(param_idx{j},["D:\Maria\OF_Excel\",".xlsx"]),'.png'])
    close all
    
end



%%
% %%%%%%%%%%%%% Ratemaps

for i =1:length(params.transcoords)

    if isempty(params.rateMap{i})
        continue
    end

map=params.rateMap{i};

figure;
PerfectCircRateMap(map,1);
 title(erase(param_idx{i},["D:\Maria\OF_Excel\",".xlsx","_"]))

print(gcf,'-dpng', '-r300',['D:\Maria\Figures\',erase(param_idx{i},["D:\Maria\OF_Excel\",".xlsx"]),'heatmaps.png'])
    close all

end


%%
%%Plot cue and circle
figure; plot(cos(0:2*pi/1000:2*pi)*101,sin(0:2*pi/1000:2*pi)*101,'-k','LineWidth',2);
 c=[params.transcue{19}(4,:);params.transcue{19}(2,:);params.transcue{19}(1,:);params.transcue{19}(3,:)];

  cue = polyshape(c);
  hold on; cue=plot(cue)
  cue.FaceColor=[0 0 0];
  cue.FaceAlpha=1;

hold on; plot(params.transcue{19}(:,1),params.transcue{19}(:,2))
axis image
axis off

%%Plot Circle Without Cue (red Star)

figure; plot(cos(0:2*pi/1000:2*pi)*101,sin(0:2*pi/1000:2*pi)*101,'-k','LineWidth',2);
  hold on; plot(cueX,cueY,'*r', 'MarkerSize',25,'LineWidth',3) 
axis image
axis off

%%

%Create example of dwell times with cue added

[zones] = createZones([0,0], 202,'fig',1,'numquad',16);
hold on; 
plot(cueX,cueY,'*r', 'MarkerSize',25,'LineWidth',3) 
hold on; 
plot(cueX*-1,cueY*-1,'*b', 'MarkerSize',25,'LineWidth',3) 

%%DWELL TIME
dwellQuad_tg1=vertcat(params.dwellQuad{contains(param_idx,'lgOF') & contains(param_idx,'Tg') & contains(param_idx,'day1')});
dwellQuad_wt1=vertcat(params.dwellQuad{contains(param_idx,'lgOF') & contains(param_idx,'Wt') & contains(param_idx,'day1')});

dwellQuad_tg2=vertcat(params.dwellQuad{contains(param_idx,'lgOF') & contains(param_idx,'Tg') & contains(param_idx,'day2')});
dwellQuad_wt2=vertcat(params.dwellQuad{contains(param_idx,'lgOF') & contains(param_idx,'Wt') & contains(param_idx,'day2')});

%Histograms for quadrant dwell time
bincenter=[0:360/16:360]-22.5;
bincenter=bincenter(2:end);

%Histograms for quadrant dwell time
f=figure; f.Color=[1 1 1];
subplot(1,2,1)
bar(sum(dwellQuad_wt1,1)/sum(dwellQuad_wt1(:)),'FaceColor','k','FaceAlpha',.5);
% title('Day 1 "Cue" Quadrant Dwell Time for TG (red) and WT (gray)')
ylim([0 .6])
xlabel('Quadrant')
ylabel('Normalized Dwell Time')
 hold on 
 bar(sum(dwellQuad_tg1,1)/sum(dwellQuad_tg1(:)),'FaceColor','r','FaceAlpha',.5);
 set(gca,'FontSize',24,'FontWeight','bold','FontName','Helvetica','LineWidth',2,'XTick',1:2:16)%'XTick',1:16,'XTickLabel',round(bincenter),
box off 

subplot(1,2,2)
bar(sum(dwellQuad_wt2,1)/sum(dwellQuad_wt2(:)),'FaceColor','k','FaceAlpha',.5);
% title('Day 2 "No Cue" Quadrant Dwell Time for TG (red) and WT (gray)')
ylim([0 .6])
xlabel('Quadrant')
ylabel('Normalized Dwell Time')
 hold on 
 bar(sum(dwellQuad_tg2,1)/sum(dwellQuad_tg2(:)),'FaceColor','r','FaceAlpha',.5);
  set(gca,'FontSize',24,'FontWeight','bold','FontName','Helvetica','LineWidth',2,'XTick',1:2:16); %'XTick',1:16,'XTickLabel',round(bincenter),
box off 

%%
























%%%%%%%%%%%%%%CODE GRAVEYARD%%%%%%%%%%%%%%%%%%
% % % % % %%"BEE SWARM" PLOT of Path Length
% % % % %
% % % % % param_idx=params.PathName;
% % % % % addpath('F:\Users\BClarkLab\GoogleDrive\MatlabDir\BClarkToolbox\Analysis\Visualize\plotSpread')
% % % % %
% % % % % %%Creates variables of each group and each day
% % % % %     pathL_tg1=vertcat(params.pathL{contains(param_idx,'lgOF') & contains(param_idx,'Tg') & contains(param_idx,'day1')});
% % % % %     pathL_wt1=vertcat(params.pathL{contains(param_idx,'lgOF') & contains(param_idx,'Wt') & contains(param_idx,'day1')});
% % % % %
% % % % %     pathL_tg2=vertcat(params.pathL{contains(param_idx,'lgOF') & contains(param_idx,'Tg') & contains(param_idx,'day2')});
% % % % %     pathL_wt2=vertcat(params.pathL{contains(param_idx,'lgOF') & contains(param_idx,'Wt') & contains(param_idx,'day2')});
% % % % %
% % % % % %%Varargins for the function and to customize the graph
% % % % %         %%Important Note: shownMW creates a line for either mean median or both, '1' creates for
% % % % %         %%both(used in line 116)
% % % % % data=[pathL_tg1; pathL_wt1;pathL_tg2;pathL_wt2];
% % % % % group=[ones(1,12) (ones(1,12)+1)]';
% % % % % catIdx=[group; group];
% % % % % dayIdx = ['Day1';'Day2']; %%labels x axis
% % % % % xValues_idx=[1;2];%the values of the x axis
% % % % % distBeesIdx=[zeros(24,1); ones(24,1)];
% % % % % yLabel_idx= ['Path Length Score'];
% % % % %
% % % % %
% % % % % %%Plotting
% % % % % h=figure;
% % % % % h=plotSpread(data,'categoryIdx',catIdx,'xNames',dayIdx,'xValues',xValues_idx,...
% % % % %     'categoryColors',{'g',[.5,.5,.5]},'DistributionIdx',distBeesIdx,'yLabel',yLabel_idx,'showMM',0,...
% % % % %     'distributionMarkers','.');
% % % % % set(h,'MarkerSize',25)
% % % % % title('Path Length')
% % % % % legend('Tg','Wt')
% % % % %
% % % % % box off
% % % % %


