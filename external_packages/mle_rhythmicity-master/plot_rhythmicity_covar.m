function [h] = plot_rhythmicity_covar( everything, plot_axis, varargin )
%PLOT_RHYTHMICITY_COVAR Plots the results of rhythmicity_covar
%   Plots the results of rhythmicity_covar. Lags can only be plotted
%   against one covariate at a time. Called by rhythmicity_covar to plot.
%   It plots three subplots (a scatter plot of lags versus the covariate,
%   a density plot of lags versus the covariate, and the mle distribution
%   of lags) within the current axis.
%
%   plot_rhythmicity_covar(everything, plot_axit);
%   h = plot_rhythmicity_covar(everything, plot_axit);
%   plot_rhythmicity_covar(everything, plot_axit, ...);
%   h = plot_rhythmicity_covar(everything, plot_axit, ...);
%
% Example: Plotting first and second covariates
% 
%   [ params, confidence, stats, everything ] = rhythmicity_covar( data,
%   session_duration, covars, 'plotit', false);
%   subplot(1,2,1);
%   plot_rhythmicity_covar(everything, 1);
%   subplot(1,2,2);
%   plot_rhythmicity_covar(everything, 2);
%
% INPUT
%   everything - everything output from rhythmicity_covar
%   plot_axis - The covariate to plot against
%
% PARAMETERS
%   CURVERES (25) - the number of points in the curves on the crest plots
%   LAGBINS (60) - The number of bins across lags
%   COVARBINS (25) - The number of bins across the covariate
%   smth (1) - The standard deviation of the smoothing kernel (bins)
%   covar_label ('') - Label for the covariate axis
%
% RETURNS
%   h - The axis handles for each of the subplots.
%
% See also cif_generator, epoch_data, rhythmicity_pdf, rhythmicity_covar
% 
% Copyright 2015-2016 Trustees of Boston University
% All rights reserved.
%
% This file is part of mle_rhythmicity revision 2.0. The last committed
% version of the previous revision is the SHA starting with 93862ac...
%
% This code has been freely distributed by the authors under the BSD
% license (http://opensource.org/licenses/BSD2-Clause). If used or
% modified, we would appreciate if you cited our papers:
%
% Climer JR, DiTullio R, Newman EL, Hasselmo ME, Eden UT. (2014),
% Examination of rhythmicity of extracellularly recorded neurons in the
% entorhinal cortex. Hippocampus, 25:460-473. doi: 10.1002/hipo.22383.
%
% Hinman et al., Multiple Running Speed Signals in Medial Entorhinal
% Cortex, Neuron (2016). http://dx.doi.org/10.1016/j.neuron.2016.06.027

% Parse inputs
ip = inputParser;
ip.addParamValue('CURVERES',25);
ip.addParamValue('LAGBINS',60);
ip.addParamValue('COVARBINS',25);
ip.addParamValue('smth',1);
ip.addParamValue('covar_label','');
ip.parse(varargin{:});
for i = fields(ip.Results)'
    eval([i{1} '=ip.Results.' i{1} ';']);
end

if isnumeric(covar_label)
   covar_label = everything.covar_labels{covar_label}; 
end

PARAMS = {'r','tau','b','c','f','s','v'};% List of parameters - used for easier coding
% If noskip, drop s from PARAMS
if everything.noskip
    PARAMS = PARAMS(~ismember(PARAMS,'s'));
end
if ~everything.refractory_rise
    PARAMS = PARAMS(~ismember(PARAMS,'v'));
end

% Subplotting
ax = gca;
pos = get(gca,'position');
HOLD = ishold;
if ~HOLD
    cla;
    axis off
end

% Scatterplot of lags versus the covariate
h(1)=subplot('position',pos.*[1 1 1 9/34]+[0 25/34*pos(4) 0 0]);
plot(everything.lags_list,everything.covars_list(:,plot_axis+1),'ok','markerfacecolor','k','markersize',1);
set(gca,'XLim',[0 everything.max_lag],'YLim',round(quantile(everything.covars_list(:,plot_axis+1),[0.01 0.95])));
xlabel('Lag (s)');ylabel(covar_label);
title('Scatter');

% Crests math
if ~isequal(everything.cov.f,0)&&ismember(plot_axis,everything.cov.f)
A = ones(CURVERES,1)*mean(everything.covars_list(:,everything.cov.f+1));
A(:,everything.cov.f==plot_axis) = linspace(...
    quantile(everything.covars_list(:,plot_axis+1),0.01)...
    ,quantile(everything.covars_list(:,plot_axis+1),0.95)...
    ,CURVERES);
f=A*everything.phat(sum(cellfun(@(x)numel(everything.cov.(x)),PARAMS(1:4)))+(1:numel(everything.cov.f)))';
A = A(:,everything.cov.f==plot_axis);
else % constant frequency
  A = linspace(...
    quantile(everything.covars_list(:,plot_axis+1),0.01)...
    ,quantile(everything.covars_list(:,plot_axis+1),0.95)...
    ,CURVERES);
  f = repmat(everything.phat(sum(cellfun(@(x)numel(everything.cov.(x)),PARAMS(1:3)))+1),...
      size(A));
end

% Plot the crests
hold on;
arrayfun(@(i)plot(i./f,A,'r--','LineWidth',2),0:ceil(max(f)*everything.max_lag));
hold off;

% Binning for density plot
COVARBINS = linspace(...
    quantile(everything.covars_list(:,plot_axis+1),0.01)...
    ,quantile(everything.covars_list(:,plot_axis+1),0.95)...
    ,COVARBINS+1);
LAGBINS = linspace(0,everything.max_lag,LAGBINS+1);

B = histcn([everything.lags_list everything.covars_list(:,plot_axis+1)],...
    LAGBINS,COVARBINS)';% Counts of lags in each bin
B = B(1:end-1,1:end-1);
B = diag(1./sum(B,2))*B;% Normalize so the sum is 1 in each row
B=imfilter(B, fspecial('gaussian',[5 5], smth), 'replicate');% Smooth

h(2)=subplot('position',pos.*[1 1 1 9/34]+[0 25/68*pos(4) 0 0]);
imagesc(LAGBINS,COVARBINS,B);
set(gca,'YDir','normal');

% Plot crests
hold on;
arrayfun(@(i)plot(i./f,A,'w--','LineWidth',2),0:ceil(max(f)*everything.max_lag));
hold off;
xlabel('Lag (s)');ylabel(covar_label);
title('Lag density');

% Relative probability from model
h(3)=subplot('position',pos.*[1 1 1 9/34]);

% 2X dense bins for this plot
COVARBINS = linspace(...
    quantile(everything.covars_list(:,plot_axis+1),0.01)...
    ,quantile(everything.covars_list(:,plot_axis+1),0.95)...
    ,(numel(COVARBINS)-1)*2);
LAGBINS = linspace(0,everything.max_lag,(numel(LAGBINS)-1)*2);
[x,y] = meshgrid(LAGBINS,COVARBINS);

% Get cif
% Generate CIFs and CIF integral functions
if everything.noskip
    if everything.refractory_rise
        [cif_fun, ~] = cif_generator('noskip_rise');
    else
    [cif_fun, ~] = cif_generator('noskip');
    end
else
    if everything.refractory_rise
         [cif_fun, ~] = cif_generator('full_rise');
    else
    [cif_fun, ~] = cif_generator('full');
    end
end
% Make design matrix
A = ones(numel(x),1)*mean(everything.covars_list);
A(:,plot_axis+1) = y(:);

% Plot
temp = reshape(...
    passall(@(varargin)...
    cif_fun(x(:),varargin{:})...
    ,covar_wrapper(everything.phat,everything.cov,A,true))...
    ,size(x));
temp = diag(sum(temp,2).^-1)*temp;
imagesc(LAGBINS,COVARBINS,temp);
set(gca,'YDir','normal');
xlabel('Lag (s)');ylabel(covar_label);
title('Relative likelihood');

end

