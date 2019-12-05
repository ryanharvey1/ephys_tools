% identify_multimodal_cells is intended to automatically classify and
% measure cells that fire maximally at multiple head angles (within yaw
% plane). 

%Code is in the process of being developed. 

%To be addressed: 

% 1. find and prevent reason for fitvmdist failures and create permanent
% patch. 

% 2. Kuiper test is too sensitive? Find more reliable test of circular
% uniformity. 

% 3. Find beter method of  hyperparameter optimization for fitvmdist, 
% perhaps instead, create new approach using methods described in Olson et al., 2016. 

% 4. Create cross-validation method to verify model selection. 

%set parameters
numDist = 4;
bin_centers = 3:6:360-3;

temp = load('F:\ClarkP30_Recordings\ProcessedData\LB07_S20190429174041.mat');


 for cell = 1:length(temp.Spikes)
    
    for ii=1:3
        
        %Check spike quality 
        
        num_spikes = size(temp.Spikes{cell},1);
        tuning = temp.hdTuning{cell,ii};
        peak_FR = max(tuning);
        
        if num_spikes < 100 || peak_FR < 2 || ii>4
            final_model(cell,ii) = NaN;
            continue
        end
        
        
        % get tuning curve to create samples
        tuning = temp.hdTuning{cell,ii};
        
        dir_prop = round(tuning*10); % use number to create tuning curve
        
        samples_deg=[]; 
        for d = 1:size(bin_centers,2)
            samples_deg = [samples_deg; repmat(bin_centers(1,d),dir_prop(1,d),1)];
        end
        
        samples=wrapToPi(deg2rad(samples_deg));
        
        % create uniform circular distribution 
        vmm = VonMisesMixture(1, 0, 0);  % Initialize model.
        uniform = vmm.random(size(samples,1));          % Draw n random samples to match data.
        
        % Check if data is uniformally distributed
         [pval, kpval, Kpval(ii,1)] = circ_kuipertest(uniform, samples, 100, 1);
        
        
        if pval < .05 %&& pval > .05
            fig=figure;
            fig.Color=[1 1 1]
            
            % evaluate models with 1 through 4 components
            for iii=1:numDist
                
                try 
                fittedVmm= fitmvmdist(samples, iii, ...
                    'MaxIter', 250); % Set maximum number of EM iterations to 300
                
                catch 
                    final_model(cell,ii) = NaN;
                    break
                end
                fit_measures(iii)=fittedVmm.logLikelihood;
            
                angles = linspace(-pi, pi, size(samples,1))';
                likelihoodsFitted = fittedVmm.pdf(angles);
                
                subplot(1,numDist,iii)
                postprocessFigures.plot_HD_tuning(temp,ii,cell)
                hold on;
                Polarplot=polar(angles, rescale(likelihoodsFitted,0,max(temp.hdTuning{cell,ii})));
                set(Polarplot,'linewidth',2,'color','r');
                title([num2str(fit_measures(iii)),' numcomponents=  ',num2str(fittedVmm.nComponents)])
                
                % initialize final_model
                final_model(cell,ii) = 0;
            end
            
            if isnan( final_model(cell,ii))
                break
            end
            
            winner = find([diff(fit_measures),0]./abs(fit_measures)>=.2)+1;
            if isempty(winner)
                winner = 1;
            end
            final_model(cell,ii) = winner(1);

        end

    end
 end
    