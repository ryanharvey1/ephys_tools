% identify_multimodal_cells is intended to classify and
% measure cells that fire maximally at multiple head angles (within
% horizontal plane).

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
nfolds = 5;
%Test data
data = load('F:\ClarkP30_Recordings\ProcessedData\LB07_S20190429174041.mat');


% for cell = 1:length(temp.Spikes) % loop across cells
%     
%     for ii = 1:3 % loop across sessions
%         
        %Check spike quality
        cell = 31; ii = 1;
        
        num_spikes = size(data.Spikes{cell},1);
        tuning = data.hdTuning{cell,ii};
        peak_FR = max(tuning);
        
        % not including low firing cells
        if num_spikes < 100 || peak_FR < 2 || ii>4
            final_model(cell,ii) = NaN;
%             continue
        end
        
        % get tuning curve to create samples
        tuning = data.hdTuning{cell,ii};
        
        dir_prop = round(tuning*10); % use number to create tuning curve
        
        samples_deg=[];
        for d = 1:size(bin_centers,2)
            samples_deg = [samples_deg; repmat(bin_centers(1,d),dir_prop(1,d),1)];
        end
        
        samples = wrapToPi(deg2rad(samples_deg));
        
        % create uniform circular distribution
        vmm = VonMisesMixture(1, 0, 0);  % Initialize model.
        uniform = vmm.random(size(samples,1));          % Draw n random samples to match data.
        
        % Create incidies for 10-fold cross validation
        indices = crossvalind('Kfold',samples,nfolds);
        
        
        for i = 1 : nfolds
            test_samples = samples(indices == i,1);
            train_samples = samples(indices ~= i,1);
           
            disp(['test set ',num2str(i),' train with all other sets'])

            res_train{i,1} = runit(train_samples, numDist);
            res_test{i,1} = runit(test_samples, numDist);
            
        end
        
%         if isnan( final_model(cell,ii))
% %             break
%         end
%         
%     end
% end


function res = runit(samples, maxOrder)

    %% model fit
    % (may take a while)
    for model_num = 1:maxOrder
        do = 1;
        while do
            try
                GMModel(model_num).res = fitmvmdist(samples,model_num);
                do = 0;
            catch
                do = 1;
            end
        end
    end
    
    %% Calculate BIC
    BIC = arrayfun(@(x) log(length(samples))*x.res.nComponents - 2*x.res.logLikelihood, GMModel);

    
     %% get parameters
     for i = 1:maxOrder
         mu{i} = arrayfun(@(x) x.mu, GMModel(i).res,'UniformOutput',0);
         k{i} =  arrayfun(@(x) x.componentProportion, GMModel(i).res,'UniformOutput',0);
     end
     
     %% package
     res.mu = mu;
     res.BIC = BIC;
     res.k = k;
     res.fits = GMModel;
     
end

% %% BIC plot
%     if ifPlot == 1
%         figure();hold on
%         plot(BIC)
%         xlabel('Model Order')
%         ylabel('BIC')
%         [minVal, minArg] = min(BIC);
%         plot(minArg,minVal,'ro','MarkerFaceColor','r')
% 
%     
%         %% Plot model fits
%         figure
%         for i = 1:5
%             axesHandle(i) = subplot(2,5,i);
%             
%             % tricsky subplotting polar plots
%             polarAxesHandle(i) = polaraxes('Units',axesHandle(i).Units,'Position',axesHandle(i).Position);
%             delete(axesHandle(i));
%             
%             % back to normal
%             binEdges = linspace(-pi,pi,30);
%             counts = histcounts(data,binEdges);
%             counts = counts / sum(counts);
%             polarhistogram(polarAxesHandle(i), 'BinEdges',binEdges, 'BinCounts', counts)
% 
%             hold on
%             evalPts = -pi:0.01:pi; evalPts = evalPts(:);
%             ft = pdf(GMModel(i).res,evalPts);
%             ft = ft / 4;  % Honestly no clue why this works ... 
%             polarplot(evalPts,ft,'LineWidth',3);
%             title(['Model Order = ' num2str(i)])
%         end
%     end