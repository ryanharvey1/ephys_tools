%% Description
% This will plot the results of all the preceding analyses: the model
% performance, the model-derived tuning curves, and the firing rate tuning
% curves.

%% plot the tuning curves

% create x-axis vectors
hd_vector = 2*pi/n_dir_bins/2:2*pi/n_dir_bins:2*pi - 2*pi/n_dir_bins/2;
ego_vector = ego_bins; 
theta_vector = hd_vector;
% speed_vector = 2.5:50/n_speed_bins:47.5;
% dist_vector = dist_bins;
% plot the tuning curves

figure(1)
subplot(3,4,1)
ratemap = pos_curve;
imAlpha=ones(size(ratemap));
imAlpha(isnan(ratemap))=0;
imagesc(ratemap,'AlphaData',imAlpha);
colormap(viridis(255))
colorbar
axis xy; axis off; box off; axis image;
title('Position')

subplot(3,4,2)
plot(hd_vector,hd_curve,'k','linewidth',3)
box off
axis([0 2*pi -inf inf])
xlabel('direction angle (rad)')
title('Head direction')

% subplot(3,4,3)
% plot(ego_vector,ego_curve,'k','linewidth',3)
% box off
% axis([0 360 -inf inf])
% xlabel('Egocentric Bearning (deg)')
% title('Egocentric Bearing')

subplot(3,4,3)
plot(theta_vector,theta_curve,'k','linewidth',3)
xlabel('Phase of Theta')
axis([0 2*pi -inf inf])
box off
title('Theta')

%% compute and plot the model-derived response profiles

% show parameters from the full model
param_full_model = param{1};

% pull out the parameter values
pos_param = param_full_model(1:n_pos_bins^2);
hd_param = param_full_model(n_pos_bins^2+1:n_pos_bins^2+n_dir_bins);
ego_param = param_full_model(n_pos_bins^2+n_dir_bins+1:n_pos_bins^2+n_dir_bins+n_ego_bins);
% dist_param = param_full_model(numel(param_full_model)-n_ego_bins+1:numel(param_full_model));

% compute the scale factors
% NOTE: technically, to compute the precise scale factor, the expectation
% of each parameter should be calculated, not the mean.
scale_factor_pos = mean(exp(ego_param))*mean(exp(hd_param))*mean(exp(dist_param))*50;
scale_factor_hd = mean(exp(ego_param))*mean(exp(pos_param))*mean(exp(dist_param))*50;
scale_factor_ego = mean(exp(pos_param))*mean(exp(hd_param))*mean(exp(dist_param))*50;
scale_factor_dist = mean(exp(ego_param))*mean(exp(hd_param))*mean(exp(pos_param))*50;

% compute the model-derived response profiles
pos_response = scale_factor_pos*exp(pos_param);
hd_response = scale_factor_hd*exp(hd_param);
% ego_response = scale_factor_ego*exp(ego_param);
% dist_response = scale_factor_dist*exp(dist_param);

% plot the model-derived response profiles
subplot(3,4,5)
imagesc(reshape(pos_response,20,20)); 
colorbar
axis xy; axis off; box off; axis image;
subplot(3,4,6)
plot(hd_vector,hd_response,'k','linewidth',3)
axis([0 2*pi -inf inf])
xlabel('direction angle')
box off
subplot(3,4,7)
plot(ego_vector,ego_response,'k','linewidth',3)
axis([0 360 -inf inf])
xlabel('Egocentric Bearing')
box off
subplot(3,4,8)
plot(dist_vector,dist_response,'k','linewidth',3)
xlabel('Distance from Boundary')
axis([0 dist_range(1,2)*2 -inf inf])
box off

% make the axes match
subplot(3,4,1)
caxis([min(min(pos_response),min(pos_curve(:))) max(max(pos_response),max(pos_curve(:)))])
subplot(3,4,5)
caxis([min(min(pos_response),min(pos_curve(:))) max(max(pos_response),max(pos_curve(:)))])

subplot(3,4,2)
axis([0 2*pi min(min(hd_response),min(hd_curve)) max(max(hd_response),max(hd_curve))])
subplot(3,4,6)
axis([0 2*pi min(min(hd_response),min(hd_curve)) max(max(hd_response),max(hd_curve))])

subplot(3,4,3)
axis([0 360 min(min(ego_response),min(ego_curve)) max(max(ego_response),max(ego_curve))])
subplot(3,4,7)
axis([0 360 min(min(ego_response),min(ego_curve)) max(max(ego_response),max(ego_curve))])

subplot(3,4,4)
axis([0 dist_range(1,2)*2 min(min(dist_response),min(dist_curve)) max(max(dist_response),max(dist_curve))])
subplot(3,4,8)
axis([0 dist_range(1,2)*2 min(min(dist_response),min(dist_curve)) max(max(dist_response),max(dist_curve))])


%% compute and plot the model performances

% ordering:
% pos&hd&spd&theta / pos&hd&spd / pos&hd&th / pos&spd&th / hd&spd&th / pos&hd /
% pos&spd / pos&th/ hd&spd / hd&theta / spd&theta / pos / hd / speed/ theta

LLH_increase_mean = mean(LLH_values);
LLH_increase_sem = std(LLH_values)/sqrt(numFolds);
if ~isnan(selected_model)
figure(1)
subplot(3,4,9:12)
errorbar(LLH_increase_mean,LLH_increase_sem,'ok','linewidth',3)
hold on
plot(selected_model,LLH_increase_mean(selected_model),'.r','markersize',25)
plot(0.5:15.5,zeros(16,1),'--b','linewidth',2)
hold off
box off
set(gca,'fontsize',20)
set(gca,'XLim',[0 16]); set(gca,'XTick',1:15)
set(gca,'XTickLabel',{'PHED','PHE','PHD','PED','HET','PH','PE','PD','HE',...
    'HD','ED','P','H','E','D'});
legend('Model performance','Selected model','Baseline')
   
end
