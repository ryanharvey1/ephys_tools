% inspect_decoding_results

files=dir('D:\Projects\PAE_PlaceCell\decoding\*.mat');

control={'RH13','RH14','LS21','LS23','LE2821','LE2823','LEM3116','LEM3120'};
pae={'RH11','RH16','LS17','LS19','LE2813','LEM3124'};

for i =1:length(files)
    load(fullfile(files(i).folder,files(i).name))
    if exist('mse','var')
        id{i,:}=files(i).name;
        all_mse(i,:)=mse;
        all_R2(i,:)=R2;
        all_mean_error(i,:)=mean_error;
        all_med_error(i,:)=med_error;
        clear mse R2 mean_error med_error
    end
end
g1_idx=contains(id,control);
g2_idx=contains(id,pae);
g1_decode=[nanmean(all_mse(g1_idx),2),nanmean(all_R2(g1_idx),2),...
    nanmean(all_mean_error(g1_idx),2),nanmean(all_med_error(g1_idx),2)];
g2_decode=[nanmean(all_mse(g2_idx),2),nanmean(all_R2(g2_idx),2),...
    nanmean(all_mean_error(g2_idx),2),nanmean(all_med_error(g2_idx),2)];

figure
AllStatsca1=stat_plot(g1_decode,g2_decode,{'Sacc','PAE'},{'mse','R2','mean error','median error'},'plots',2)
darkBackground(gcf,[0.2 0.2 0.2],[0.7 0.7 0.7])

figure
subplot(4,1,1)
errorbar(1:size(all_mse,1),nanmean(all_mse,2),nanstd(all_mse')*sqrt(1/10+1/(10-1)),'o');hold on
scatter(1:size(all_mse,1),nanmean(all_mse,2),20,sum([g1_idx,g2_idx*2],2))
ylabel('mse')
grid on

subplot(4,1,2)
plot([1,size(all_R2,1)],[0,0],'r');hold on
errorbar(1:size(all_R2,1),nanmean(all_R2,2),nanstd(all_R2')*sqrt(1/10+1/(10-1)),'o');hold on
scatter(1:size(all_R2,1),nanmean(all_R2,2),20,sum([g1_idx,g2_idx*2],2))
ylabel('R2')
grid on

subplot(4,1,3)
errorbar(1:size(all_mean_error,1),nanmean(all_mean_error,2),nanstd(all_mean_error')*sqrt(1/10+1/(10-1)),'o');hold on
scatter(1:size(all_mean_error,1),nanmean(all_mean_error,2),20,sum([g1_idx,g2_idx*2],2))
ylabel('mean error')
grid on

subplot(4,1,4)
errorbar(1:size(all_med_error,1),nanmean(all_med_error,2),nanstd(all_med_error')*sqrt(1/10+1/(10-1)),'o');hold on
scatter(1:size(all_med_error,1),nanmean(all_med_error,2),20,sum([g1_idx,g2_idx*2],2))
ylabel('med error')
grid on
xlabel('sessions')

colormap cool
darkBackground(gcf,[0.2 0.2 0.2],[0.7 0.7 0.7])

% influence of number of cells
for i=1:length(id)
    load(erase(id{i},'_decode_results'),'Spikes')
    ncells(i)=length(Spikes);
end
figure
subplot(2,2,1)
scatter(ncells,nanmean(all_mse,2),20,sum([g1_idx,g2_idx*2],2))
lsline
xlabel('number of cells')
ylabel('mse')
grid on

subplot(2,2,2)
scatter(ncells,nanmean(all_R2,2),20,sum([g1_idx,g2_idx*2],2))
lsline
xlabel('number of cells')
ylabel('R2')
grid on

subplot(2,2,3)
scatter(ncells,nanmean(all_mean_error,2),20,sum([g1_idx,g2_idx*2],2))
lsline
xlabel('number of cells')
ylabel('mean error')
grid on

subplot(2,2,4)
scatter(ncells,nanmean(all_med_error,2),20,sum([g1_idx,g2_idx*2],2))
lsline
xlabel('number of cells')
ylabel('med error')
grid on

colormap cool
darkBackground(gcf,[0.2 0.2 0.2],[0.7 0.7 0.7])


% influence of session duration
for i=1:length(id)
    load(erase(id{i},'_decode_results'),'session_duration')
    sess_duration(i)=session_duration(1);
end
figure
subplot(2,2,1)
scatter(sess_duration,nanmean(all_mse,2),20,sum([g1_idx,g2_idx*2],2))
lsline
xlabel('session duration') 
ylabel('mse')
grid on

subplot(2,2,2)
scatter(sess_duration,nanmean(all_R2,2),20,sum([g1_idx,g2_idx*2],2))
lsline
xlabel('session duration') 
ylabel('R2')
grid on

subplot(2,2,3)
scatter(sess_duration,nanmean(all_mean_error,2),20,sum([g1_idx,g2_idx*2],2))
lsline
xlabel('session duration') 
ylabel('mean error')
grid on

subplot(2,2,4)
scatter(sess_duration,nanmean(all_med_error,2),20,sum([g1_idx,g2_idx*2],2))
lsline
xlabel('session duration') 
ylabel('med error')
grid on

colormap cool
darkBackground(gcf,[0.2 0.2 0.2],[0.7 0.7 0.7])

% group=sum([g1_idx,g2_idx*2],2);
% mse=nanmean(all_mse,2);
% [h,atab,ctab,stats] = aoctool(mse,ncells,group);
% [h,atab,ctab,stats] = aoctool(mse,sess_duration,group)
% 
% [r,p]=corr(nanmean(all_mse,2),ncells')
% 
% tbl = table(nanmean(all_mse,2),sess_duration',ncells',...
%     sum([g1_idx,g2_idx*2],2),...
%     'VariableNames',{'mse','sess_duration','ncells','group'});
% tbl.group=categorical(tbl.group);
% gzlm=fitglm(ds,'Response ~ Group+Covariate+Group*Covariate','Distribution','poisson');
% 
% fit = fitglm(tbl,'mse ~ group*sess_duration+group*ncells')
% anova(fit)