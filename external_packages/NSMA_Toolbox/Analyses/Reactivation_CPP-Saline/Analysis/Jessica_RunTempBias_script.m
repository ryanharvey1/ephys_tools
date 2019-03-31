% RunTempBias script    Run it from the Analysis dir as pwd (present working dir) dir!
% 
% PL May 2004 for datsets obeying the session object syntax 
% This version is adapted to Jessica's 2004 CPP data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Select Session Group %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% select group of sessions for current run
group = 'CPP'            % sessions are selected (filtered) by this session.group parameter
%group = 'Saline'         % sessions are selected (filtered) by this session.group parameter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Set RUN NAME and OUTPUT DIR  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

outPrefix = [group '_'];
outputDir = ['Out\TempBias_30Spikes_Binsize200ms_' group '_15min'];

% create diary and make sure output dir exists
if ~exist(outputDir,'dir')
    mkdir(outputDir);
end
diary off;                  % in case it is on from a bombed previous run
diary([outputDir filesep outPrefix, 'run_', strrep(datestr(datenum(now),0),':','-'),'.txt']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Set RUN PARAMTERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%set limits for cuts on cell firing properties
SpikeCountLowerLimit_S1 = 30 % sleep 1 lower spike count limit
SpikeCountLowerLimit_S2 = 30 % sleep 2 lower spike count limit
SpikeCountLowerLimit_B  = 30 % behavior epoch lower spike count limit
rateLim = 2.2               % mean firing rate upper limit (Hz) for cells during behavior for inclusion in the analysis

%set restriction intervals for sleep 1 and 2 epochs
sleep1_offset_sec = 60      % 60 sec offset before end of sleep 1. 
sleep1_interval_sec = -900  % restrict sleep1 to 10 min interval before end of sleep1 minus offset
sleep2_offset_sec = 180     % give animal 3 min to start sleep2 
sleep2_interval_sec = 900   % restrict sleep2 to 10 min interval after offset

% paramters for CCG-matrix and generation
exclude_sameTT = 1          % exclude cell pairs from same tetrode  
binsize_S = 10         % sleep binsize:  5 for 100ms windows or 2.5 for 50 ms windows
binsize_B = 10         % behavior binsize: binsize for CCGs in ms
numbins = 100         % number of bins (also timelags) for CCGs
normal = 1;           % if normal = 1, correct for outliers in bias (often resulting in 1 or -1 for bias) 
RATHC_DATA = 1;       % This tells the script that the data is rat hippocampus data
                      % which alters certain options such as size of the cell arrays used

% load the session list
sessList = SessionList01;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% global variables
global RATHC_DATA
global normal
global irat
global isess
global namesize
global sessname
global ratlist

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
sessListPerRat
ratlist

sessname = {};
% create sessnames cell array
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
%fh = figure;

% initialize pool cell arrays
NCellsTotal = 0;
sess_count = 0;
cs1_pool = {};
cm_pool = {};
cs2_pool = {};

%%%%%%%%%%%%%%%%%%%%%%%%% Loop over animals %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for irat = 1:length(ratlist)
    
    NCellsPerRat = 0;    % Counter of number of used cells in analysis (after all cuts)

    % Initialize the cell array to hold the "all" data. Placed in the for
    %  loop because these cell arrays get large so they are removed after the
    %  data for each rat is completed and then re-initialized
    if RATHC_DATA
        sysCCGS1all{1} = [];
        sysCCGBall{1} = [];
        sysCCGS2all{1} = [];
    else        
        sysCCGS1all{4,4}=[];
        sysCCGBall{4,4}=[];
        sysCCGS2all{4,4}=[];
        sysULS1all{4,4}=[];
        sysULBall{4,4}=[];
        sysULS2all{4,4}=[];
    end %if  

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
        
        nPairs = SameNTrodeMatrix_CountPairs(sameTT);
        fprintf('\n RAT %s, Session %s:\n session has N = %d cells and nPairs = %d valid pairs after cuts\n\n',ses.animal,ses.name,length(SB),nPairs); 
    
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

        % count number of contributing cells and sessions
        N = [length(SB),0,0,0];  
        NCellsPerRat = NCellsPerRat + N(1);
        nsess = nsess + 1;            % count sessions contributing to the pool for this animal
        sess_count = sess_count + 1;  % count  all sessions contributing to the pool 
        
        % make CCGmats for each epoch
        CCGS1 = doCCG (SS1, binsize_S, numbins,sameTT);
        CCGB = doCCG (SB, binsize_B, numbins,sameTT);
        CCGS2 = doCCG (SS2, binsize_S, numbins,sameTT);
        
        % pass epoch S's into vectorizer to get CCGs by systemcombo (10 possible)  
        % sort out within electrode CCGs for ccgmat, and CI mats
        [sysCCGS1split{isess}, allsz1] = CCGsess2sys(N, CCGS1);
        [sysCCGBsplit{isess}, allszB]  = CCGsess2sys(N, CCGB);
        [sysCCGS2split{isess}, allsz2] = CCGsess2sys(N, CCGS2);
        
        % combine across sessions for each syscombo
        % Rest1
        for isys = 1:length(sysCCGS1split{isess}) % iterate for each syscombo 
            for jsys = 1:length(sysCCGS1split{isess}) 
                sysCCGS1all{isys,jsys} = [sysCCGS1all{isys,jsys}; sysCCGS1split{isess}{isys,jsys}];
            end % for jsys
        end % for isys
        % Behavior
        for isys = 1:length(sysCCGBsplit{isess}) % iterate for each syscombo 
            for jsys = 1:length(sysCCGBsplit{isess}) 
                sysCCGBall{isys,jsys} = [sysCCGBall{isys,jsys}; sysCCGBsplit{isess}{isys,jsys}];
            end % for jsys
        end % for isys
        % Rest 2
        for isys = 1:length(sysCCGS2split{isess}) % iterate for each syscombo 
            for jsys = 1:length(sysCCGS2split{isess}) 
                sysCCGS2all{isys,jsys} = [sysCCGS2all{isys,jsys}; sysCCGS2split{isess}{isys,jsys}];
            end % for jsys
        end % for isys
        
        
    end % for isess
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End loop over sessions
    
    if nsess == 0
        % ignore current animal
        fprintf('\n ===> RAT %s has contributed 0 sessions to pool and is skipped!!\n\n',ses.animal); 
        continue;
    end

    CCGbias  % compute the bias of all CCG variables (all and split)
    
    % CLEAN UP THE BIAS VECTORS
    CCGRemoveNaNs   % Remove all NaNs from the Bias results (usually caused by 0/0)
                    % which results from all of the bins 31:50 and 52:71 being 0
    
    if (normal)              
        CCGRemoveOnes   %Remove all of the -1 or 1 values from the Bias, this usually occurs
    end %if             % when only one or two bins on one side has any values such that, e.g.
                        % the sum of -200ms to 0ms <<<<< sum of 0ms to 200ms.
    
    % Compute the correlation between S1bias - Bbias and S2bias - Bbias 
    ABM = [S1allbias{1} Ballbias{1} S2allbias{1}];
    ABM(:,4) = (ABM(:,3) - ABM(:,1));                   % Sleep 2 - Sleep 1 
    
    % for EV 
    S1B_allbiascorr = diag(corrcoef(ABM(:,1),ABM(:,2)),1);
    S2B_allbiascorr = diag(corrcoef(ABM(:,3),ABM(:,2)),1);
    S2B_allbiassub = diag(corrcoef(ABM(:,4),ABM(:,2)),1);
    S1S2_allbiascorr = diag(corrcoef(ABM(:,1),ABM(:,3)),1);
    
    
    %Partial Correlation
    norm = sqrt((1 - S1B_allbiascorr^2) * (1 - S1S2_allbiascorr^2));    
    S2B_allbiaspart = (S2B_allbiascorr - S1B_allbiascorr*S1S2_allbiascorr)/ norm;
    % compute the reversed partial corrleation = sqrt( explained variance (EV)) with role of S1 and S2 interchanged for m,s2|s1
    norm_rev = sqrt((1 - S2B_allbiascorr^2) * (1 - S1S2_allbiascorr^2));    
    S2B_allbiaspart_rev = (S1B_allbiascorr - S2B_allbiascorr*S1S2_allbiascorr)/ norm_rev;
    
    fprintf('\n=============================================================================================\n');
    fprintf('\n RAT %s, %d sessions pooled together into %d cells:\n EV = %f, EV_rev = %f \n',...
        ses.animal, nsess, NCellsPerRat, S2B_allbiaspart, S2B_allbiaspart_rev); 
    fprintf('\n=============================================================================================\n');
   
    
    % This part builds a vector of the bias correlation results for each dataset, I'm not
    % sure that this is set to work for multiple arrays yet.
    S1B_sessbiascorr = [];
    S2B_sessbiascorr = [];
    if exist('S1sessbias','var')
		for x = 1:length(S1sessbias)
            if length(S1sessbias{x}) > 1
                S1B_sessbiascorr(x) = diag(corrcoef(S1sessbias{x},Bsessbias{x}),1);
                S2B_sessbiascorr(x) = diag(corrcoef(S2sessbias{x},Bsessbias{x}),1);
            else
                S1B_sessbiascorr(x) = NaN;
                S2B_sessbiascorr(x) = NaN;
            end
        end
    end
    
    animalSessionList = sessListPerRat{irat};
    NCellsTotal = NCellsTotal + NCellsPerRat;

    % save the data to the a file that is called "RAT#_bias_results.mat"
    filename = [outputDir filesep ratlist{irat} '_' num2str(binsize_S *20) 'ms_bias_results'];
    save(filename,'sys*','*allbia*','*sessbia*','ABM','NCellsPerRat','animalSessionList','nsess');
    
    % Clear out the memory and move on to the next rat
    clear *sessbia*; 
    clear *allbia*; 
    clear CCG*; 
    clear S1* S2* SS* SB*; 
    clear post pre; 
    clear *isnotnan; 
    clear sameTT;
    clear sys*; 
    clear tfiles*; 
    clear ts*; 
    %pack;
       
    
end% for irat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End loop over rats



% save figure
%filename = [outputDir filesep outPrefix 'BiasHistograms_AllAnimals_AllSessions' num2str(Btime_binsize) '_' limname 'Hz_Results'];
%saveas(fh,filename,'fig');

diary off

