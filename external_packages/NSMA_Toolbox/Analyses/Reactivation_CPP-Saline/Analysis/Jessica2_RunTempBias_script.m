% RunTempBias script 2
% Removed all dependency on Kari's and Jason's code
% 
% PL April 2003 for datsets obeying the session object syntax 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Select Session Group %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% select group of sessions for current run
runNr = 12
%group = 'CPP'            % sessions are selected (filtered) by this session.group parameter
group = 'Saline'         % sessions are selected (filtered) by this session.group parameter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Set RUN NAME and OUTPUT DIR  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
outPrefix = [group '_'];
outputDir = ['T:\CPP1\JessicasData\Analysis\out\TempBias_' group '_15min'];

% create diary and make sure output dir exists
if ~exist(outputDir,'dir')
    disp(['making new directory: ' outputDir ]);
    mkdir(outputDir);
end
fclose('all');              % close all files in case one is left open from a bombed previous run
diary off;                  % in case it is on from a bombed previous run
runNrStr = sprintf('%03d',runNr);
runFileName = [outputDir, filesep, 'Run',runNrStr,'_', outPrefix, strrep(datestr(datenum(now),0),':','-')];
runFileName = strrep(runFileName,' ','_')
diary([runFileName '.txt']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Open CSV files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fpCSV1 = fopen([runFileName '_sessions.csv'],'wt');
fpCSV2 = fopen([runFileName '_rats.csv'],'wt');
% write CSV headers
fprintf(fpCSV1,['Rat,','Session,','Group,','Length of Sleep 1 (sec),', 'length of Behavior (sec),', 'length of Sleep 2 (sec),','# Cells before cuts,','# Cells after cuts,','# contributing Pairs','\n']);
fprintf(fpCSV2,['Rat,','Group,','# Sessions,','# Cells per rat,','# Pairs before cuts,', '# Pairs after NaNs and 1s removal,', '# Pairs after area cut,', 'Partial Corr,','Rev. Partial Corr','\n']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Set RUN PARAMTERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%set limits for cuts on cell firing properties
SpikeCountLowerLimit_S1 = 30 % sleep 1 lower spike count limit
SpikeCountLowerLimit_S2 = 30 % sleep 2 lower spike count limit
SpikeCountLowerLimit_B  = 30 % behavior epoch lower spike count limit
rateLim = 2.2               % mean firing rate upper limit (Hz) for cells during behavior for inclusion in the analysis

%set restriction intervals for sleep 1 and 2 epochs
sleep1_offset_sec = 0      % 60 sec offset before end of sleep 1. 
sleep1_interval_sec = -900  % restrict sleep1 to 10 min interval before end of sleep1 minus offset
sleep2_offset_sec = 0     % give animal 3 min to start sleep2 
sleep2_interval_sec = 900   % restrict sleep2 to 10 min interval after offset

% paramters for CCG-matrix and generation
exclude_sameTT = 1          % exclude cell pairs from same tetrode  
binsize_S = 10              % sleep binsize:  5 for 100ms windows or 2.5 for 50 ms windows
binsize_B = 10              % behavior binsize: binsize for CCGs in ms
numbins = 100               % number of bins (also timelags) for CCGs
normal = 1                  % if normal = 1, correct for outliers in bias (often resulting in 1 or -1 for bias) 

% cut on area of crosscorrelogram of cell pairs
area_cut = 0
area_threshold = 0         % Run 000
%area_threshold = exp(8)    % Run 001
%area_threshold = exp(10)    % Run 002

% load the session list
sessList = SessionList01;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % global variables
% global normal
% global irat
% global isess
% global namesize
% global sessname
% global ratlist

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
CCGS1all = {};
CCGBall = {};
CCGS2all = {};
CCGS1 = {};
CCGB = {};
CCGS2 = {};
areaMall = [];
clear bias* pre* post*;


%%%%%%%%%%%%%%%%%%%%%%%%% Loop over animals %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for irat = 1:length(ratlist)
    
    NCellsPerRat = 0;    % Counter of number of used cells in analysis (after all cuts)

    CCGS1rat = {};
    CCGBrat = {};
    CCGS2rat = {};

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
        NCellsPreCuts = length(tfiles_S1{isess});
        fprintf('\n RAT %s, Session %s:\n session has N = %d cells prior to cuts\n\n',ses.animal,ses.name,NCellsPreCuts); 

        
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
        N = length(SB);  
        NCellsPerRat = NCellsPerRat + N;
        nsess = nsess + 1;            % count sessions contributing to the pool for this animal
        sess_count = sess_count + 1;  % count  all sessions contributing to the pool 
        
        % make cross-correlogram cell arrays for each epoch
        [CCGS1{irat,isess}, bincenters] = MakeCrossCorrCellArray(SS1, binsize_S, numbins, sameTT);
        [CCGB{irat,isess}, bincenters]  = MakeCrossCorrCellArray(SB,  binsize_B, numbins, sameTT);
        [CCGS2{irat,isess}, bincenters] = MakeCrossCorrCellArray(SS2, binsize_S, numbins, sameTT);
        NParisPreCuts = length(CCGS1{irat,isess});
        
        % combine pairs for each rat across sessions 
        CCGS1rat = [CCGS1rat, CCGS1{irat,isess} ];
        CCGBrat = [CCGBrat, CCGB{irat,isess} ];
        CCGS2rat = [CCGS2rat, CCGS2{irat,isess} ];
        
        % save session info to CVS
        fprintf(fpCSV1,['%s,'     ,'%s,'   ,'%s,','%f,'           ,'%f,'          ,'%f,'           ,'%d,'        ,'%d,','%d'           ,'\n'],...
                        ses.animal,ses.name,group,SS1_interval_sec,SB_interval_sec,SS2_interval_sec,NCellsPreCuts, N   , NParisPreCuts );

    end % for isess
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End loop over sessions
    
    if nsess == 0
        % ignore current animal
        fprintf('\n ===> RAT %s has contributed 0 sessions to pool and is skipped!!\n\n',ses.animal); 
        continue;
    end

     % compute the bias of all CCG variables 
     norm = 1;          % compute normalized (1) or unnormalized (0) bias
     width = 19;        % nunmber of bins in pre- and post areas of cross-correlograms
    [biasS1rat, preS1rat, postS1rat] = getCCGBias(CCGS1rat, norm, width);
    [biasBrat,  preBrat,  postBrat]  = getCCGBias(CCGBrat,  norm, width);
    [biasS2rat, preS2rat, postS2rat] = getCCGBias(CCGS2rat, norm, width);
    NPairsPreCuts = length(biasS1rat);
    
    % Make a Bias-Matrix with each row: [biasS1, biasB, biasS2, biasS1-biasS1] and one row per pair of cells 
    ABM = [biasS1rat biasBrat biasS2rat];
    ABM(:,4) = (ABM(:,3) - ABM(:,1));                   % Sleep 2 - Sleep 1 
    preM = [preS1rat preBrat preS2rat];
    postM = [postS1rat postBrat postS2rat];
    
    
    % CLEAN UP THE BIAS MATRIX:
    % Remove all NaNs from the Bias matrix (usually caused by 0/0)
    % which results from all of the bins 31:50 and 52:71 being 0
	ix_rows = find(sum(isnan(ABM),2));   
	ABM(ix_rows,:) = [];
    preM(ix_rows,:) = [];
    postM(ix_rows,:) = [];
    [nr, nc] = size(ABM);
    fprintf('\n Removing %d pairs (rows) with NaNs from the Bias Matrix, leaving %d pairs.' ,length(ix_rows),nr)
    if (normal)              
        % Remove all of the -1 or 1 values from the Bias, this usually occurs
        % when only one or two bins on one side has any values such that, e.g.
        % the sum of -200ms to 0ms <<<<< sum of 0ms to 200ms.
		ix_rows = find(sum(abs(ABM(:,1:3))==1,2));
		ABM(ix_rows,:) = [];
        preM(ix_rows,:) = [];
        postM(ix_rows,:) = [];
        [nr, nc] = size(ABM);
        fprintf('\n Removing %d pairs (rows) with 1s from the Bias Matrix, leaving %d pairs.' ,length(ix_rows),nr)
    end %if             
    NPairsPostNaNsCut = nr;
    
    % AREA CUT:
    % Remove all rows from the Bias matrix which correspond to pairs with area in Crosscorrlogram < area_cut
    if area_cut
        areaM = preM+postM;
        ix_rows = find(sum(areaM < area_threshold,2));
        areaM(ix_rows,:) = [];
		ABM(ix_rows,:) = [];
        preM(ix_rows,:) = [];
        postM(ix_rows,:) = [];
        [nr, nc] = size(ABM);
        fprintf('\n Removing %d pairs (rows) below area_threshold %f from the Bias Matrix, leaving %d pairs.' ,length(ix_rows),area_threshold,nr)
    end %if
    NPairsPostAreaCut = nr;
    
    % accumulate areaM for all rats
    if area_cut
        areaMall = [areaMall; areaM];
    else
        areaMall = [];
    end
    
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
        ses.animal, nsess, NCellsPerRat, S2B_allbiaspart(1), S2B_allbiaspart_rev(1)); 
    fprintf('\n=============================================================================================\n');
   
        
    animalSessionList = sessListPerRat{irat};
    NCellsTotal = NCellsTotal + NCellsPerRat;

    % save rat info to CSV
    fprintf(fpCSV2,['%s,','%s,','%d,','%d,'        ,'%d,'         , '%d,'            , '%d,'            , '%f,'          ,'%f'                ,'\n'], ...
            ses.animal,  group, nsess, NCellsPerRat, NPairsPreCuts, NPairsPostNaNsCut, NPairsPostAreaCut, S2B_allbiaspart(1), S2B_allbiaspart_rev(1));

    % save the data to the a file that is called "RAT#_bias_results.mat"
    filename = [outputDir, filesep, 'Run',runNrStr,'_', 'Bias_results_', ratlist{irat}, '_', num2str(binsize_S *20), 'ms'];
    fprintf('\n writing variables to file %s in dir %s' ,filename,pwd)
    save(filename,'bias*','pre*','post*','CCG*','*allbia*','ABM','NCellsPerRat','animalSessionList','nsess');
           
    
end% for irat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End loop over rats


% save figure
%filename = [outputDir filesep outPrefix 'BiasHistograms_AllAnimals_AllSessions' num2str(Btime_binsize) '_' limname 'Hz_Results'];
%saveas(fh,filename,'fig');

fclose('all');
diary off



