function plot_glm_HD_res(data,session,glm_res,varargin)

p = inputParser;
addOptional(p,'n_circ_bins',18,@isnumeric)
addOptional(p,'n_speed_bins',10,@isnumeric)
parse(p,varargin{:})


% Initialize model hyperparameters
n_dir_bins = p.Results.n_circ_bins;
n_speed_bins = p.Results.n_speed_bins;


% boxSize = length (in cm) of one side of the square box
boxSize=data.maze_size_cm(session);
frames=data.frames(data.frames(:,1)>data.events(1,session) &...
    data.frames(:,1)<data.events(2,session),:); 

% convert xy to cm 
frames(:,2)=rescale(frames(:,2),0,boxSize);
frames(:,3)=rescale(frames(:,3),0,boxSize);

% Grab position and head direction 
posx_c = frames(:,2);
posy_c = frames(:,3);
direction = deg2rad(frames(:,4));

% for velocity filtering
[~,~,speed] = speed_map(posx_c,posy_c,n_speed_bins);
% take out times when the animal ran >= 50 cm/s
too_fast = find(speed >= 100);
direction(too_fast) = [];

% compute tuning curves for position, head direction, speed, and theta phase
[hd_curve] = compute_1d_tuning_curve(direction,glm_res.smooth_fr,n_dir_bins,0,2*pi);

% create x-axis vectors
hd_vector = 2*pi/n_dir_bins/2:2*pi/n_dir_bins:2*pi - 2*pi/n_dir_bins/2;

% plot the tuning curves
fig=figure;fig.Color=[1 1 1];
subplot(1,3,1)
plot(hd_vector,hd_curve,'k','linewidth',3)
box off
axis([0 2*pi -inf inf])
xlabel('direction angle')
title('Head direction')

% Get predicted tuning curves from models
hd_response = response_profile(glm_res.param);

% plot the model-derived response profiles
subplot(1,3,2)
plot(hd_vector,hd_response,'k','linewidth',3)
xlabel('direction angle')
title('Predicted Tuning')
axis([0 2*pi -inf inf])
box off


% make the axes match
subplot(1,3,1)
axis([0 2*pi min(min(hd_response),min(hd_curve)) max(max(hd_response),max(hd_curve))])
subplot(1,3,2)
axis([0 2*pi min(min(hd_response),min(hd_curve)) max(max(hd_response),max(hd_curve))])

% Loglikelihood values
LLH_values = reshape(glm_res.test(:,3),10,1);
LLH_increase_mean = mean(LLH_values);
LLH_increase_sem = std(LLH_values)/sqrt(10);

subplot(1,3,3)
errorbar(LLH_increase_mean,LLH_increase_sem,'ok','linewidth',3)
hold on
if ~isnan(glm_res.best_model)
    plot(glm_res.best_model,LLH_increase_mean(glm_res.best_model),'.r','markersize',25)
end
hold off
box off
grid on
set(gca,'fontsize',12)
set(gca,'XLim',[.5 1.5]); set(gca,'XTick',1)
set(gca,'XTickLabel',{'HD'});
legend('Model performance','Selected model','Baseline')
ylabel('Log-Likelihood (bits/spike)')

end


function [tuning_curve] = compute_1d_tuning_curve(variable,fr,numBin,minVal,maxVal)

%bin it
var_vec = linspace(minVal,maxVal,numBin+1);
tuning_curve = nan(numBin,1);

% compute mean fr for each bin
for n = 1:numBin
    tuning_curve(n) = nanmean(fr(variable >= var_vec(n) & variable < var_vec(n+1)));
    
    if n == numBin
        tuning_curve(n) = nanmean(fr(variable >= var_vec(n) & variable <= var_vec(n+1)));
    end
    
end


end

function[hd_response] = response_profile(param)

% compute the model-derived response profiles
hd_response = exp(param)*33;

end