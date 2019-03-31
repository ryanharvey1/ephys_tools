function [params, confidence, stats, everything] = mle_rhythmicity( data, session_duration, varargin )
%MLE_RHYTHMICITY Maximum likelihood estimation for rhythmicity parameters
%   Using a conditional intensity function for the distribution of lags in
%   the autocorrelogram window, the set of parameters that maximize the
%   likelihood of the data are identified.
%
%   [...] = mle_rhythmicitiy(data, session_duration, [PARAMETER], [value]...);
%
%   INPUT:
%      data: Can either be spike timestamps or a cell array of lags
%      session_duration: Duration of session
%
%   PARAMETERS
%       max_lag (0.6): Examination window (seconds)
%       plotit (true): If true, plots histogram and distribution estimates
%       noskip (false): If true, omits skipping from the model
%       refractory_rise (false): If true, adds an exponential rise time
%           constant (v)
%       f_range ([1 13]): Range of possible frequencies for the rhythmicity.
%       alpha (0.05): Cutoff for confidence intervals
%       t_axis (linspace(0,maxlag,61)): Axis for plotting
%       runs (3): Number of runs, takes the best answer of all.
%
%   OUTPUT
%       params: A struct containing the parameter estimates
%           a - Relative magnitude of the rhythmicity, ranges from 0
%               (no rhythm) to 1 (maximally rhythmic)
%           tau (log10(sec)) - Exponential falloff of the distrubtion
%           b - Baseline
%           c (log10(sec)) - Exponential falloff of the rhythmicity magnitude
%           f (Hz) - Frequency of the rhythmicity
%           s - Relative heights of the secondary peaks. Not present if
%               noskip=true
%           v - Half-rise time for the exponential refractory rise (sec).
%               Not present if refractory_rise=false.
%       confidence: Struct containing confidence intervals of the parameters
%       stats: Struct containing
%           deviance_rhythmic - Deviance for the rhythmicity test. If
%               noskip=false, this has 4 degrees of freedom, otherwise 3.
%           p_rhythmic - Significance of chi-squared test for rhythmicity
%           deviance_skipping - Deviance for the skipping test, with 1 degree
%               of freedom. Not present if noskip=true
%           p_rhythmic - Significance of chi-squared test for skipping. Not
%               present if noskip=true.
%       everything: Struct containing
%           LL - Log likelihood for the fit
%           LL_flat - Log likelihood for the non-rhythmic fit
%           se - Standard errors for the parameters
%           phat - Parameter estimate as a vector
%           lags - Cell array of all the lags
%           lags_list - list of all the lags
%           max_lag - Maximum lag
%           noskip - Whether noskip=true
%           phat_flat - Parameter estimates for the non-rhythmic fit as a
%               vector
%           session_duration - Duration of the session
%           LL_noskip - Log likelihood of the non-skipping fit. Not present if
%               noskip=true.
%           phat_noskip - Parameter estimates for the non-skipping fit as a
%               vector. Not present if noskip=true.
%           se_noskip - Standard errors for non-skipping parameters as a
%               vector. Not present if noskip=true.
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

%% Warnings
warning('off','stats:mle:IterLimit');

%% Parse inputs
alpha = 0;
ip = inputParser;
% This allows us to keep inputs if this was called by rhythmicity_covar
temp = dbstack;
if ismember('rhythmicity_covar',{temp.name})
    ip.KeepUnmatched = true;
end
ip.addParamValue('max_lag', 0.6);
ip.addParamValue('epochs',[]);
ip.addParamValue('epochmode','lead');
ip.addParamValue('plotit',true);
ip.addParamValue('noskip',false);
ip.addParamValue('refractory_rise',false);
ip.addParamValue('f_range',[1 13]);
ip.addParamValue('alpha',0.05);
ip.addParamValue('t_axis',[ ]);
ip.addParamValue('runs',3);
ip.parse(varargin{:});
% Put all parsed inputs into workspace
for j = fields(ip.Results)'
    eval([j{1} ' = ip.Results.' j{1} ';']);
end

if runs>1 % If we're doing multiple runs
    warning off stats:mlecov:NonPosDefHessian;
    
    % Reformat input for multiple runs
    j = true(size(varargin));
    if ~ismember('plotit',ip.UsingDefaults)
        j(find(cellfun(@(x)isequal(x,'plotit'),varargin))+[0 1]) = false;
    end
    if ~ismember('runs',ip.UsingDefaults)
        j(find(cellfun(@(x)isequal(x,'runs'),varargin))+[0 1]) = false;
    end
    
    % Run the first time
    [params, confidence, stats, everything] = mle_rhythmicity( data, session_duration, varargin{j}, 'runs', 1, 'plotit', false );
    
    % Resize for multiple runs
    params = repmat(params,[runs 1]);
    confidence = repmat(confidence,[runs 1]);
    stats = repmat(stats,[runs 1]);
    everything = repmat(everything,[runs 1]);
    
    % Do other runs
    for i=2:runs
        [params(i), confidence(i), stats(i), everything(i)] = mle_rhythmicity( data, session_duration, varargin{j}, 'runs', 1, 'plotit', false );
    end
    
    % Find the best fit
    [~,i] = max([everything.LL]);
    
    % Select everything for output
    params = params(i);
    confidence = confidence(i);
    stats = stats(i);
    everything = everything(i);
else % A single run
    
    % Format and epoch the data
    if ~iscell(data)
        if isempty(epochs)
            epochs = [0 session_duration];
        end
        [ lags, inds ] = epoch_data( data, epochs ,'epochmode', epochmode, 'max_lag', max_lag);
    end
    
    lags_list = cat(1,lags{:});
    n = numel(lags_list);
    
    %%
    % This abstract function takes in two handles for a CIF and its definite
    % integral and the parameters for the CIF, and returns the log-likelihood
    % of the set of data
    LL_fun = @(cif_fun,cif_int,varargin)infbnd(sum(log(cif_fun(lags_list,varargin{:}))-log(cif_int(max_lag,varargin{:}))));
    
    if ~refractory_rise
    [cif_fun, cif_int] = cif_generator('flat');% Non-rhythmic cif, see cif_generator
    
    [phat_flat, ci_flat] = mle(1,'logpdf',@(~,varargin)LL_fun(cif_fun,cif_int,varargin{:})...
        ,'start',[-log10(log(1/0.95)/max_lag) 0.1]...% Tau starts to complete 95% of decay by the maximum lag
        ,'lowerbound',[-inf 0]....
        ,'upperbound',[10 1]);
    LL_flat = passall(@(varargin)LL_fun(cif_fun,cif_int,varargin{:}),phat_flat);
    else
        [cif_fun, cif_int] = cif_generator('flat_rise');% Non-rhythmic cif, see cif_generator
    
    [phat_flat, ci_flat] = mle(1,'logpdf',@(~,varargin)LL_fun(cif_fun,cif_int,varargin{:})...
        ,'start',[-log10(log(1/0.95)/max_lag) 0.1 0.005]...% Tau starts to complete 95% of decay by the maximum lag
        ,'lowerbound',[-inf 0 0]....
        ,'upperbound',[10 1 max_lag/4]);
    LL_flat = passall(@(varargin)LL_fun(cif_fun,cif_int,varargin{:}),phat_flat);
    end
    % Non-skipping rhythm
    
    % First, we use a pure sinusoid to find the best frequency
    [cif_fun, cif_int] = cif_generator('pure');% Pure sinusoid, see cif_generator
    f = linspace(f_range(1),f_range(2),200);% To start, check many frequencies across available range
    [~,j] = max(arrayfun(@(f)LL_fun(cif_fun, cif_int, f),f));% Find the best one
    f = f(j);
    
    % Use builtin optimizers to converge to final best frequency
    f = mle(1,'logpdf',@(~,varargin)LL_fun(cif_fun,cif_int,varargin{:})...
        ,'start',f...
        ,'lowerbound',f_range(1)....
        ,'upperbound',f_range(2));
    
    %% Lets try to start b at the termination of the flat fit, and fix the
    % frequency above
    PopulationSize = 25;
    % tau, b, c, r
    InitialPopulation = [...
        unifrnd(-2,1,[PopulationSize,1])...tau
        unifrnd(0,0.4,[PopulationSize,1]) ...b
        unifrnd(-2,1,[PopulationSize,1])...c
        unifrnd(max(f*0.9,f_range(1)),min(f*1.1,f_range(2)),[PopulationSize,1])...f
        unifrnd(0.25,0.75,[PopulationSize,1]) ...r
        ];
    
    if refractory_rise
        InitialPopulation = [InitialPopulation(:,1:end-1)...
            unifrnd(0,max_lag/4,[PopulationSize 1])...v
            InitialPopulation(:,end)];
        [cif_fun, cif_int] = cif_generator('noskip_rise');
        lb = 0;
        ub = max_lag/4;
        
    else
        [cif_fun, cif_int] = cif_generator('noskip');
        lb = [];
        ub = [];
    end
    % Particle swarm fit
    phat_noskip = pso(...
        @(phat)-passall(@(varargin)LL_fun(cif_fun,cif_int,varargin{:}),phat)...
        ,5+refractory_rise ...
        ,[],[],[],[]...
        ,[-inf 0 -inf f_range(1) lb 0]...
        ,[10 0.5 10 f_range(2) ub 1]...
        ,[] ...
        ,psooptimset('Display','off','Generations',100,'InitialPopulation',InitialPopulation,'PopulationSize',PopulationSize,'ConstrBoundary','Reflect'...
        ...,'PlotFcns',{@psoplotswarm}...
        )...
        );
    
    % Final convergence with mle
    phat_noskip = mle(1,'logpdf',@(~,varargin)LL_fun(cif_fun,cif_int,varargin{:})...
        ,'start',phat_noskip...
        ,'lowerbound',[-inf 0 -inf f_range(1) lb 0]...
        ,'upperbound',[10 1 10 f_range(2) ub 1]);
    
    %%
    % % Fit using non-skipping distribution
    
    % Convert to A & calculate CIs, log-likelihood
    %     keyboard
    phat_noskip = [(1-phat_noskip(2))*phat_noskip(end) phat_noskip(1:end-1)];
    se_noskip = (diag(mlecov(phat_noskip,1,'logpdf',@(~,a,varargin)LL_fun(cif_fun,cif_int,varargin{:},a/(1-varargin{2})))))';
    ci_noskip = [-1;1]*se_noskip*norminv(1-alpha/2)+repmat(phat_noskip,[2 1]);
    LL_noskip = passall(@(a,varargin)LL_fun(cif_fun,cif_int,varargin{:},a/(1-varargin{2})),phat_noskip);
    
    if ~noskip % Add skipping
        
        if refractory_rise
            phat = [phat_noskip(2:end-1) 0.05 phat_noskip(end) phat_noskip(1)/(1-phat_noskip(3))];
            
            lb = 0;
            ub = max_lag/4;
            
            [cif_fun, cif_int] = cif_generator('full_rise');
        else
            phat = [phat_noskip(2:end) 0.05 phat_noskip(1)/(1-phat_noskip(3))];
            [cif_fun, cif_int] = cif_generator('full');
            
            lb = [];
            ub = [];
            
        end
        
        phat = mle(1,'logpdf',@(~,varargin)LL_fun(cif_fun,cif_int,varargin{:})...
            ,'start',phat...
            ,'lowerbound',[-inf 0 -inf f_range(1) 0 lb 0]...
            ,'upperbound',[10 1 10 f_range(2) 1 ub 1]);
        
        LL = passall(@(varargin)LL_fun(cif_fun,cif_int,varargin{:}),phat);
        
        if phat_noskip(5)*2<=f_range(2)% If doubling the frequency is still in the frequency range
            %  Fit using a doubled frequency and high skipping
            phat_half = mle(1,'logpdf',@(~,varargin)LL_fun(cif_fun,cif_int,varargin{:})...
                ,'start',[phat_noskip(2:4) phat_noskip(5)*2 0.95 phat_noskip(1)/(1-phat_noskip(3))]...
                ,'lowerbound',[-inf 0 -inf f_range(1) 0 lb 0]...
                ,'upperbound',[10 1 10 f_range(2) 1 ub 1]);
            LL_half = passall(@(varargin)LL_fun(cif_fun,cif_int,varargin{:}),phat_half);
            if LL_half>LL % If doubling the frequency is a better fit, replace
                phat = phat_half;
            end
        end
        
        % Convert to A and calculate CIs, log-likelihood
        phat = [(1-phat(2))*phat(end) phat(1:end-1)];
        se = (diag(mlecov(phat,1,'logpdf',@(~,a,varargin)LL_fun(cif_fun,cif_int,varargin{:},a/(1-varargin{2})))))';
        ci = [-1;1]*se*norminv(1-alpha/2)+repmat(phat,[2 1]);
        LL = passall(@(a,varargin)LL_fun(cif_fun,cif_int,varargin{:},a/(1-varargin{2})),phat);
        
        % Test for significant skipping
        stats.deviance_skipping = 2*(LL-LL_noskip);
        stats.p_skipping = 1-chi2cdf(stats.deviance_skipping,1);
        
        % Test for significant rhythm
        stats.deviance_rhythmic = 2*(LL-LL_flat);
        stats.p_rhythmic = 1-chi2cdf(stats.deviance_rhythmic,4);
    else % noskip=true
        % Use noskip fits
        phat = phat_noskip;
        se = se_noskip;
        ci = ci_noskip;
        LL = LL_noskip;
        
        % Test for significant rhythm
        stats.deviance_rhythmic = 2*(LL-LL_flat);
        stats.p_rhythmic = 1-chi2cdf(stats.deviance_rhythmic,3);
    end
    %% Package output
    PARAM = {'a','tau','b','c','f'};
    if ~noskip, PARAM = [PARAM {'s'}]; end
    if refractory_rise
        PARAM = [PARAM {'v'}];
    end
    for j=1:numel(PARAM)
        eval(sprintf('params.%s=phat(%i);',PARAM{j},j));
        eval(sprintf('confidence.%s=ci(:,%i);',PARAM{j},j));
    end
    
    PARAM = {'inds','epochs','LL','LL_flat','se','phat','lags','lags_list','max_lag','noskip','phat_flat','session_duration'};
    if ~noskip, PARAM = [PARAM {'LL_noskip','phat_noskip','se_noskip'}]; end
    for j=1:numel(PARAM)
        eval(sprintf('everything.%s=%s;',PARAM{j},PARAM{j}));
    end
end

%% plotting
if ~ishold
    cla
end
if plotit(1)
    if isempty(t_axis)
        t_axis = linspace(0,max_lag,61);
    else
        if ~all(diff(t_axis)-(t_axis(2)-t_axis(1))<0.0001)
            warning('Custom t_axis must be monotonically increasing. Using default 61 bins')
            t_axis = linspace(0,max_lag,61);
        end
    end
    tn = length(t_axis);
    b=histc(everything.lags_list,t_axis);
    bar(mean([t_axis(1:end-1);t_axis(2:end)]),b(1:end-1),1);
    xlim([0 max_lag]);
    hold on;
    if ~refractory_rise
        if noskip
            [cif_fun, cif_int] = cif_generator('noskip');
        else
            [cif_fun, cif_int] = cif_generator('full');
        end
    else
        if noskip
            [cif_fun, cif_int] = cif_generator('noskip_rise');
        else
            [cif_fun, cif_int] = cif_generator('full_rise');
        end
    end
    ps = passall(@(varargin)cif_fun(linspace(0,max_lag,200),varargin{2:end},varargin{1}/(1-varargin{3})),everything.phat)/...
        passall(@(varargin)cif_int(max_lag,varargin{2:end},varargin{1}/(1-varargin{3})),everything.phat);
    
   
    
    plot(linspace(0,max_lag,200),ps*max_lag*numel(everything.lags_list)/tn,'r','LineWidth',2);
    if ~refractory_rise
    [cif_fun, cif_int] = cif_generator('flat');
    else
        [cif_fun, cif_int] = cif_generator('flat_rise');
    end
    
    ps = passall(@(varargin)cif_fun(linspace(0,max_lag,200),varargin{:}),everything.phat_flat)/...
        passall(@(varargin)cif_int(max_lag,varargin{:}),everything.phat_flat);
        
    plot(linspace(0,max_lag,200),ps*max_lag*numel(everything.lags_list)/tn,'c--','LineWidth',2);
    ttl = ['\hat{a}' sprintf('=%2.2g, p',everything.phat(1)) '_{rhyth}=' sprintf('%2.2g',stats.p_rhythmic)];
    lgnd = {'data','MLE','flat'};
    if ~noskip
        if ~refractory_rise
            [cif_fun, cif_int] = cif_generator('noskip');
        else
            [cif_fun, cif_int] = cif_generator('noskip_rise');
        end
        ps = passall(@(varargin)cif_fun(linspace(0,max_lag,200),varargin{2:end},varargin{1}/(1-varargin{3})),everything.phat_noskip)/...
            passall(@(varargin)cif_int(max_lag,varargin{2:end},varargin{1}/(1-varargin{3})),everything.phat_noskip);
              
        plot(linspace(0,max_lag,200),ps*max_lag*numel(everything.lags_list)/tn,'g--');
        
        ttl = [ttl ', s=' sprintf('%2.2g',everything.phat(end)) ', p_{skip}=' sprintf('%2.2g',stats.p_skipping)];
        lgnd = [lgnd {'no-skip'}];
    end
    hold off
    legend(lgnd{:});
    %clc
    title(['$$' ttl '$$'],'Interpreter','latex');
    xlabel('Lag (s)');ylabel('Count');
end

end

