contains(param_idx,'lgOF') & contains(param_idx,'Tg');

dwellQuad_tg1=vertcat(params.dwellQuad{contains(param_idx,'lgOF') & contains(param_idx,'Tg') & contains(param_idx,'day1')});
dwellQuad_wt1=vertcat(params.dwellQuad{contains(param_idx,'lgOF') & contains(param_idx,'Wt') & contains(param_idx,'day1')});

dwellQuad_tg2=vertcat(params.dwellQuad{contains(param_idx,'lgOF') & contains(param_idx,'Tg') & contains(param_idx,'day2')});
dwellQuad_wt2=vertcat(params.dwellQuad{contains(param_idx,'lgOF') & contains(param_idx,'Wt') & contains(param_idx,'day2')});

%Histograms for quadrant dwell time
bincenter=[0:360/16:360]-22.5;
bincenter=bincenter(2:end);
%%
f=figure; f.Color=[1 1 1];
subplot(2,2,1)
bar(sum(dwellQuad_tg1,1)/sum(dwellQuad_tg1(:)),'FaceColor','r');
title('tg day 1')
ylim([0 .6])
xlabel('Quadrant Angle')
ylabel('Normalized Dwell Time')
set(gca,'XTick',[1:16],'XTickLabel',bincenter)

box off 

subplot(2,2,3)
bar(sum(dwellQuad_wt1,1)/sum(dwellQuad_wt1(:)),'FaceColor',[.5 .5 .5]);
title('wt day 1')
ylim([0 .6])
xlabel('Quadrant Angle')
ylabel('Normalized Dwell Time')
set(gca,'XTick',[1:16],'XTickLabel',bincenter)

box off 

subplot(2,2,2)
bar(sum(dwellQuad_tg2,1)/sum(dwellQuad_tg2(:)),'FaceColor','r');
title('tg day 2')
ylim([0 .6])
xlabel('Quadrant Angle')
ylabel('Normalized Dwell Time')
set(gca,'XTick',[1:16],'XTickLabel',bincenter)

box off 

subplot(2,2,4)
bar(sum(dwellQuad_wt2,1)/sum(dwellQuad_wt2(:)),'FaceColor',[.5 .5 .5]);
title('wt day 2')
ylim([0 .6])
xlabel('Quadrant Angle')
ylabel('Normalized Dwell Time')
set(gca,'XTick',[1:16],'XTickLabel',bincenter)

box off 
%%


figure; plot(cos(0:2*pi/1000:2*pi)*101,sin(0:2*pi/1000:2*pi)*101,'-k');

for i=1:length(y)
hold on; plot(x{:,i},y{:,i},'--r')
end

hold on; plot(params.transcoords{j}(:,1),params.transcoords{j}(:,2),'.b')
axis image

% hold on; comet(params.transcoords{j}(:,1),params.transcoords{j}(:,2))


%% Path Length 

pathL_tg1=vertcat(params.pathL{contains(param_idx,'lgOF') & contains(param_idx,'Tg') & contains(param_idx,'day1')});
pathL_wt1=vertcat(params.pathL{contains(param_idx,'lgOF') & contains(param_idx,'Wt') & contains(param_idx,'day1')});

pathL_tg2=vertcat(params.pathL{contains(param_idx,'lgOF') & contains(param_idx,'Tg') & contains(param_idx,'day2')});
pathL_wt2=vertcat(params.pathL{contains(param_idx,'lgOF') & contains(param_idx,'Wt') & contains(param_idx,'day2')});


