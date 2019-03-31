% RunEV script
% 
% PL Feb 2003 for datsets obeying the session object syntax 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Set RUN NAME and OUTPUT DIR  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
outPrefix = 'Saline01_'
outputDir = 'Out\EV\EV_Saline_familiar_100msecBins_FirstS2Interval'
%outPrefix = 'CPP01_'
%outputDir = 'Out\EV\EV_CPP_familiar_100msecBins_FirstS2Interval'

% create diary and make sure output dir exists
if ~exist(outputDir,'dir')
    mkdir(outputDir);
end
diary off;                  % in case it is on from a bombed previous run
diary([outputDir filesep outPrefix, 'run_', strrep(datestr(datenum(now),0),':','-'),'.txt']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Set RUN PARAMTERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% select group of sessions for current run
%group = 'CPP'            % sessions are selected (filtered) by this session.group parameter
group = 'Saline'            % sessions are selected (filtered) by this session.group parameter

%set limits for cuts on cell firing properties
SpikeCountLowerLimit_S1 = 30 % sleep 1 lower spike count limit
SpikeCountLowerLimit_S2 = 30 % sleep 2 lower spike count limit
SpikeCountLowerLimit_B  = 30 % behavior epoch lower spike count limit
rateLim = 2.2               % mean firing rate upper limit (Hz) for cells during behavior for inclusion in the analysis

%set restriction intervals for sleep 1 and 2 epochs
sleep1_offset_sec = 0      % 60 sec offset before end of sleep 1. 
sleep1_interval_sec = -600  % restrict sleep1 to 10 min interval before end of sleep1 minus offset
sleep2_offset_sec = 0       % give animal 3 min to start sleep2 
sleep2_interval_sec = 600   % restrict sleep2 to 10 min interval after offset

% paramters for R-matrix and generation
Stime_binsize = 100         % binsize of sleep Q matrices in msec
Btime_binsize = 100         % binsize of  behavior Q matrix in msec
rescale_timebins = 0        % turn on(1)/off(0) the rescaling of timebins for S1 and S2 matrices
exclude_sameTT = 1          % exclude cell pairs from same tetrode  

% load the session list
sessList = SessionList02;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% make sure cell arrays in workspace are cleared for this run
clear avgRate*;

%prepare the list of rats and sessions
sessGroup = SessionListFilter(sessList,'group',group);
ratlist = {};
ratdir = {};
sessListPerRat = {};
for i=1:length(sessGroup)
    res = strcmpi(ratlist,sessGroup{i}.animal);
    if ~sum(res)
        ratlist{end+1} = sessGroup{i}.animal;
        ratdir{end+1} = sessGroup{i}.tfiledir;
        sessListPerRat{end+1} = SessionListFilter(sessList,'group',group,'animal',sessGroup{i}.animal);
    end
end
ratlist

% create sessnames cell array
sessname = {};
for irat = 1:length(ratlist)
    for isess = 1:length(sessListPerRat{irat})
        ses = sessListPerRat{irat}{isess};
        sessname{irat,isess} = [ses.animal,'_',ses.name];
    end%for
end%for irat
sessname


% Get the dimensions of the cell array containing all the datasets 
[nRats, nSessMax] = size(sessname);
limname = strrep(num2str(rateLim),'.','_')  % needed for output filenames - do not change        

% initialize summary plots:
fh = figure;

% initialize pool cell arrays
NCellsTotal = 0;
sess_count = 0;
cs1_pool = {};
cm_pool = {};
cs2_pool = {};

%%%%%%%%%%%%%%%%%%%%%%%%% Loop over animals %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for irat = 1:length(ratlist)
    
    NCellsPerRat = 0;    % Counter of number of used cells in analysis (after all cuts)

    %For the number of columns in sessname (largest possible number of datasets)  
    working = {};
    for z = 1:nSessMax
        if (~isempty(sessname{irat,z}))     % if this cell isn't empty - it has a dataset name in it
            working{z} = sessname{irat,z};  % write the dataset name into the temp variable working
        end %if
    end
    % working
    % Now working is the list of datasets to be analyzed for the current rat
    
    clear tfiles*;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Loop over sessions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % For the number of datasets listed for this rat
    nsess = 0;
    for isess = 1:length(working)
 
        ses = sessListPerRat{irat}{isess};
        fprintf('\n\n\n======================== RAT: %s, Session: %s, Group: %s ==============\n', ses.animal,ses.name, ses.group); 

        tfiles_S1{isess} =  getTFileList(ses,'sleep1');
        tfiles_B{isess}  =  getTFileList(ses,'maze1');
        tfiles_S2{isess} =  getTFileList(ses,'sleep2');
        fprintf('\n RAT %s, Session %s:\n session has N = %d cells prior to cuts\n\n',ses.animal,ses.name,length(tfiles_S1{isess})); 

        
        SS1 = getSpikeMatrix(ses,'sleep1');
        SB  = getSpikeMatrix(ses,'maze1');
        SS2 = getSpikeMatrix(ses,'sleep2');
        
        sameTT = getSameNTrodeMatrix(ses);
        
        % restrict the S-matrices to the desired sub-intervals:
        SS1 = RestrictTheSMatrix2Interval(SS1,sleep1_offset_sec,sleep1_interval_sec);
        SS2 = RestrictTheSMatrix2Interval(SS2,sleep2_offset_sec,sleep2_interval_sec);
        
        % check intervals
        [ts0,ts1] = SessionStartEndTS(SS1);
        SS1_interval_sec = (ts1-ts0)/10000;
        fprintf('SS1: ts0 = %d, ts1 = %d, interval(sec) = %f\n',ts0,ts1,SS1_interval_sec);
        [ts0,ts1] = SessionStartEndTS(SB);
        SB_interval_sec = (ts1-ts0)/10000;
        fprintf('SB: ts0 = %d, ts1 = %d, interval(sec) = %f\n',ts0,ts1,SB_interval_sec);
        [ts0,ts1] = SessionStartEndTS(SS2);
        SS2_interval_sec = (ts1-ts0)/10000;
        fprintf('SS2: ts0 = %d, ts1 = %d, interval(sec) = %f\n\n',ts0,ts1,SS2_interval_sec);
        
        % check avg firing rates
        fprintf('Average Firing Rates prior to cuts:\n');
        avg_rate = SMatrix_getAvgFiringRates(SS1);
        fprintf('SS1: avgRate = %f +- std = %f (Hz)\n',mean(avg_rate),std(avg_rate));
        avg_rate = SMatrix_getAvgFiringRates(SB);
        fprintf('SB:  avgRate = %f +- std = %f (Hz)\n',mean(avg_rate),std(avg_rate));
        avg_rate = SMatrix_getAvgFiringRates(SS2);
        fprintf('SS2: avgRate = %f +- std = %f (Hz)\n\n',mean(avg_rate),std(avg_rate));

        % apply cuts on tfiles: compile a list of cells to be discarded for this analysis
        disp('applying cuts on tfiles:');
        badindx = [];        
        for i = 1:length(SS1)
            secs = (EndTime(SB{i})-StartTime(SB{i}))/10000;

            if (length(data(SS1{i})) < SpikeCountLowerLimit_S1)     % must be at least 30 spikes during S1
                badindx = [ badindx ; i ];
                fprintf('The tfile %s had less than %d spikes during S1 \n',tfiles_S1{isess}{i},SpikeCountLowerLimit_S1);
                
            elseif (length(data(SS2{i})) < SpikeCountLowerLimit_S2)     % must be at least 30 spikes during S2
                badindx = [ badindx ; i ];
                fprintf('The tfile %s had less than %d spikes during S2 \n',tfiles_S2{isess}{i},SpikeCountLowerLimit_S2);
                
            elseif (length(data(SB{i})) < SpikeCountLowerLimit_B)        %must be at least 90 spikes during behavior
                badindx = [ badindx ; i];
                fprintf('The tfile %s had less than %d spikes during Behavior \n',tfiles_B{isess}{i},SpikeCountLowerLimit_B);
            
            elseif ((length(data(SB{i}))/secs) > rateLim)  % Check mean firing rate during behavior
                 badindx = [ badindx ; i];
                 fprintf('The tfile %s had a mean rate greater than 2.2 Hz during behavior \n',tfiles_B{isess}{i});
            end %if
        end
        
        % discard the list of cells not passing above criteria
        if ~isempty(badindx)
            disp('Found tfiles that did not pass all cuts')
            fprintf('%d cells will be removed from the analysis.\n',length(badindx));           
            tfiles_S1{isess}(badindx) = [];
            tfiles_S2{isess}(badindx) = [];
            tfiles_B{isess}(badindx) = [];
            SS1(badindx) = [];
            SB(badindx) = [];
            SS2(badindx) = [];
            sameTT(:,badindx) = [];
            sameTT(badindx,:) = [];
        end%if
        
        
        N = length(SB);  
        NCellsPerRat = NCellsPerRat + N;  
        nPairs = SameNTrodeMatrix_CountPairs(sameTT);
        fprintf('\n RAT %s, Session %s:\n session has N = %d cells and nPairs = %d valid pairs after cuts\n\n',ses.animal,ses.name,N,nPairs); 
    
        if nPairs == 0
            % ignore current session
            fprintf('\n ===> RAT %s, Session %s:\n this session has not enough pairs (Pairs = %d) and is skipped!!\n\n',ses.animal,ses.name,nPairs); 
            continue;
        end
        
        % check avg firing rates again after cuts
        fprintf('Average Firing Rates after cuts:\n');
        avg_rate = SMatrix_getAvgFiringRates(SS1);
        avgRateSS1{irat,isess} = [mean(avg_rate),std(avg_rate)];
        fprintf('SS1: avgRate = %f +- std = %f (Hz)\n',mean(avg_rate),std(avg_rate));
        avg_rate = SMatrix_getAvgFiringRates(SB);
        avgRateSB{irat,isess} = [mean(avg_rate),std(avg_rate)];
        fprintf('SB:  avgRate = %f +- std = %f (Hz)\n',mean(avg_rate),std(avg_rate));
        avg_rate = SMatrix_getAvgFiringRates(SS2);
        avgRateSS2{irat,isess} = [mean(avg_rate),std(avg_rate)];
        fprintf('SS2: avgRate = %f +- std = %f (Hz)\n\n',mean(avg_rate),std(avg_rate));


        % compute correlation vectors (flattended R-matrices) for current session: SS1, SS2 and SB and sameTT
        [cs1,cm,cs2] = ...
            RMatrixCorrCoeffs(SS1,SB,SS2,sameTT,Stime_binsize,Btime_binsize,exclude_sameTT,rescale_timebins);

        nsess = nsess + 1;            % count sessions contributing to the pool for this animal
        sess_count = sess_count + 1;  % count  all sessions contributing to the pool 
        cs1_pool{end+1} = cs1;
        cm_pool{end+1} = cm;
        cs2_pool{end+1} = cs2;
       
        
        % calculate Explained variance
        [ev,ev_rev] = ExplainedVariance(cs1,cm,cs2);
        fprintf('\n RAT %s, Session %s:\n EV = %f, EV_rev = %f \n\n',ses.animal,ses.name,ev,ev_rev); 

        
        % store results per session
        [s1corrs_hist, bbs1] = hist(cs1,[-.1:.01:.5]);
        [mcorrs_hist, bbm] = hist(cm,[-.1:.01:.5]);
        [s2corrs_hist, bbs2] = hist(cs2,[-.1:.01:.5]);
        % make summary figure of corrs
        figure(fh);
        subplot(nRats, nSessMax, (irat-1)*nSessMax + isess);
        hold on;
        plot(bbs1,s1corrs_hist,'b');
        plot(bbm, mcorrs_hist,'r');
        plot(bbs2,s2corrs_hist,'g');
        axis tight;
        title (strrep(sessname{irat,isess},'_','-'));
        hold off;
        drawnow;

		filename = [outputDir filesep outPrefix sessname{irat,isess}  '_EV_' num2str(Btime_binsize) '_' limname 'Hz_Results'];
		corrfilename = [outputDir filesep outPrefix sessname{irat,isess} '_correlation_' num2str(Btime_binsize) '_' limname 'Hz_lists'];
		save(corrfilename,'cs1', 'cm', 'cs2');
		save(filename, 'r*', '*corrs_hist', 'ev', 'ev_rev', 'bb*', 'tfiles*','*_interval_sec','nPairs');
        
        %store results for all rats/sessions
        EV_results{irat,isess} = {sessname{irat,isess},ev,ev_rev};
                
        
    end % for isess
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End loop over sessions
    
    if nsess == 0
        % ignore current animal
        fprintf('\n ===> RAT %s has contributed 0 sessions to pool and is skipped!!\n\n',ses.animal); 
        continue;
    end

    
    % compute EVs of the all the sessions of the current animal pooled together:
    cs1 = cat(1,cs1_pool{end-nsess+1:end});
    cm  = cat(1,cm_pool{end-nsess+1:end});
    cs2 = cat(1,cs2_pool{end-nsess+1:end});
    [ev, ev_rev] = ExplainedVariance(cs1,cm,cs2);
    fprintf('\n=============================================================================================\n');
    fprintf('\n RAT %s, %d sessions pooled together into %d cells and %d pairs:\n EV = %f, EV_rev = %f \n',...
        ses.animal, nsess, NCellsPerRat, length(cs1), ev, ev_rev); 
    fprintf('\n=============================================================================================\n');

    EV_animal_results = EV_results(end,1:nsess);
    animalSessionList = sessListPerRat{irat};
    filename = [outputDir filesep outPrefix 'EV_' ratlist{irat} '_AllSessions' num2str(Btime_binsize) '_' limname 'Hz_Results'];
    save(filename,'animalSessionList','EV_animal_results','NCellsPerRat','ev','ev_rev','nsess','cs1','cm','cs2')
    
    NCellsTotal = NCellsTotal + NCellsPerRat;
    
end% for irat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End loop over rats


% compute EVs of the all sessions of all animals pooled together:
cs1 = cat(1,cs1_pool{:});
cm  = cat(1,cm_pool{:});
cs2 = cat(1,cs2_pool{:});
[ev, ev_rev] = ExplainedVariance(cs1,cm,cs2);
fprintf('\n=============================================================================================\n');
fprintf('=============================================================================================\n');
fprintf('\n ALL RATS, %d sessions pooled together into %d cells and %d pairs:\n EV = %f, EV_rev = %f \n',...
    sess_count, NCellsTotal, length(cs1), ev, ev_rev); 
fprintf('\n=============================================================================================\n');
fprintf('=============================================================================================\n');


filename = [outputDir filesep outPrefix 'EV_AllAnimals_AllSessions' num2str(Btime_binsize) '_' limname 'Hz_Results'];
save(filename,'sessListPerRat','avgRate*','EV_results','NCellsTotal','ev','ev_rev','nsess','cs1','cm','cs2')

% save figure
filename = [outputDir filesep outPrefix 'CorrHistograms_AllAnimals_AllSessions' num2str(Btime_binsize) '_' limname 'Hz_Results'];
saveas(fh,filename,'fig');

diary off

