function [num_components , mean_prom, score] = get_nComponents(data,session, cell, fig)

% Input: 
%       data: data structure specified in ephys tools. Must contain
%       subfield hdTuning. 
%           note. hdTuning is a 1 x 60 tuning curve reflecting binned spike
%           count over 6 degree bins normalized by occupancy. Each row in
%           data.hdTuning reflects the calculated tuning curve for a given
%           cell. 
%       session: scalar indicating session number (acts as index to select
%       column in data.hdTuning. 
%
%       cell: scalar indicating the cell index (acts as index to select row
%       in data.hdTuning. 
%
%       fig: boolean indicating,  0 = no figure, 1 = figure. Figure shows tuning curve, smoothed
%       tuning curve for which the autocorr is produced from and the peaks
%       obtained from the autcorr. 

% Output: 
%       num_components: number of peaks identified in circular autocor 
%
%       mean_prom: average height of the peaks 
%
%       score: peak high minus average of troughs


% Code adapted from hdbidrectionalityscore (https://github.com/kevin-allen)

% create several example distributions to test below function (bimodal,
% unimodal, trimodal using mvmdist functions)

% % bimodal von mises
% p = [0.5; 0.5];         % Mixture weights.
% mu = [-pi/2; pi/2];     % Component means.
% kappa = [5; 1];        % Concentration parameters of components.
%
% vmm = VonMisesMixture(p, mu, kappa);
% bimodal = vmm.random(10000);          % Draw 10000 random samples
%
% % unimodal von mises
% p = [1];         % Mixture weights.
% mu = [-pi/2];     % Component means.
% kappa = [3];        % Concentration parameters of components.
%
% vmm = VonMisesMixture(p, mu, kappa);
% unimodal = vmm.random(10000);          % % Draw 10000 random samples.
%
% % trimodal von mise
% p = [1/3; 1/3; 1/3];         % Mixture weights.
% mu = [-pi/2; 0; pi/2];     % Component means.
% kappa = [10; 5; 2];        % Concentration parameters of components.
%
% vmm = VonMisesMixture(p, mu, kappa);
% trimodal = vmm.random(10000);          % % Draw 10000 random samples
%
%
% bin_centers = -pi:pi/30:pi
% hdTuning_test(1,:)=histcounts(unimodal,bin_centers);
% hdTuning_test(2,:)=histcounts(bimodal,bin_centers);
% hdTuning_test(3,:)=histcounts(trimodal,bin_centers);

%%

% smooth tuning curve

hdTuning = data.hdTuning{cell, session};

smooth_tune = smoothdata([hdTuning hdTuning hdTuning],'gaussian',6); %smooth over 36 degrees
smooth_tune = smooth_tune(1, 61:120);


x = 1:60; %for plotting

% Let's get a circular autocorr to determine when peaks occur relative
% to each other.
ac = cconv(smooth_tune,conj(fliplr(smooth_tune)),size(smooth_tune,2)); % get a circular autocorrelation

% find peaks within the autocor starting at the center
[peak, prominence] = islocalmax(ac,'MaxNumExtrema',5,'MinSeparation',6,...
    'MinProminence',1); % set max extrema to 5 and min separation to 9 (6 *6 = 36 degrees)


% find peaks that are not at beginning or end
locs = x(peak);
pks = ac(peak);
prom_pks = prominence(locs);

locs(locs == 1 | locs == 60) = [];
pks(locs == 1 | locs == 60) = [];


if isempty(locs)
    
    num_components = 1;
    mean_prom = NaN;
    score = 0;
    
    if fig
        % Plot the above
        figure;
        subplot(2,2,1)
        plot(hdTuning,'-k','LineWidth',2)
        title(['raw hdTuning: ',num2str(data.spikesID.TetrodeNum{cell}),num2str(data.spikesID.CellNum(cell))])
        subplot(2,2,2)
        plot(smooth_tune,'-r','LineWidth',2)
        title('smooth hdTuning')
        subplot(2,2,3.5)
        plot(x,ac,'LineWidth',2); hold on
        plot(x,ac,x(peak),ac(peak),'r*')
        title(['circular autocorr: num components: ', num2str(num_components),'score: ',num2str(unique(score))])
        
    end
    return
else
    
    num_components = size(locs,2) + 1; % bimodal = 1 peak, trimodal = 2 peaks,
    %                                   so add 1 to number of identified peaks
    %                                   to find number of component estimates.
    mean_prom = nanmean(prominence);
    
end

% 
% Below is an attempt to replicate the code from Kevin Allen's github.
% Only works with bimodal cells.

% find highest peak

[max_peak, ~]= max (pks); % highest peak
peak_loc = x(ac == max_peak); % highest peak location

% find troughs
[trough,~] = islocalmin(ac,'MaxNumExtrema',5,'MinSeparation',6);


% find troughs before and after peak
locs_trough = x(trough);

if size(peak_loc,2) == 1 % only one peak
    trough_before = ac(locs_trough(locs_trough < peak_loc));
    trough_after  = ac(locs_trough(locs_trough > peak_loc));
    
    % get the lowest trough before and after the trough
    
    min_before = min(trough_before);
    min_after  = min(trough_after);
    
    score = (max_peak - mean([min_before; min_after]));
    
else
    
    for k = 1:size(peak_loc,2)
        trough_before = ac(locs_trough(locs_trough < peak_loc(k)));
        trough_after  = ac(locs_trough(locs_trough > peak_loc(k)));
        
        % get the lowest trough before and after the trough
        
        min_before = min(trough_before);
        min_after  = min(trough_after);
        
        score (k,1) = (max_peak - mean([min_before; min_after]));
    end
    
end

if fig
    % Plot the above
    figure;
    subplot(2,2,1)
    plot(hdTuning,'-k','LineWidth',2)
    title(['raw hdTuning: ',num2str(data.spikesID.TetrodeNum{cell}),num2str(data.spikesID.CellNum(cell))])
    subplot(2,2,2)
    plot(smooth_tune,'-r','LineWidth',2)
    title('smooth hdTuning')
    subplot(2,2,3.5)
    plot(x,ac,'LineWidth',2); hold on
    plot(x,ac,x(peak),ac(peak),'r*')
    plot(x,ac,x(trough),ac(trough),'b*')
    title(['circular autocorr: num components: ', num2str(num_components),'score: ',num2str(unique(score))])
    
end

end
