% RunTempBias Parameter Dependence script 
%
% Study dependence of CPP/Saline Temporal Bias Reactivation (slopes of S1-B/S2-B regressions) 
% in dependence of AREA cut paramter.
%
% April 22, 2003: Removed all dependency on Kari's and Jason's code
% 
% PL April 2003 for datsets obeying the session object syntax 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Select Session Group %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% select group of sessions for current run
runNr = 30
group = 'CPP'            % sessions are selected (filtered) by this session.group parameter
%group = 'Saline'         % sessions are selected (filtered) by this session.group parameter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Set RUN NAME and OUTPUT DIR  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
outPrefix = [group '_AreaCutDependence_'];
outputDir = ['D:\CPP1\JasonCPPSleepStudy\jason\Out\TempBias_15min_NoRateCut'];

% create diary and make sure output dir exists
if ~exist(outputDir,'dir')
    mkdir(outputDir);
end
fclose('all');              % close all files in case one is left open from a bombed previous run
diary off;                  % in case it is on from a bombed previous run
runNrStr = sprintf('%03d',runNr);
runFileName = [outputDir, filesep, 'Run',runNrStr,'_', outPrefix, strrep(datestr(datenum(now),0),':','-')]
diary([runFileName '.txt']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Open CSV files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fpCSV1 = fopen([runFileName '_sessions.csv'],'wt');
fpCSV2 = fopen([runFileName '_rats.csv'],'wt');
fpCSV3 = fopen([runFileName '_runs.csv'],'wt');
% write CSV headers
fprintf(fpCSV1,['Run,','Rat,','Session,','Group,','Length of Sleep 1 (sec),', 'length of Behavior (sec),', 'length of Sleep 2 (sec),','# Cells before cuts,','# Cells after cuts,','# contributing Pairs','\n']);
fprintf(fpCSV2,['Run,','Rat,','Group,','# Sessions,','# Cells per rat,','# Pairs before cuts,', '# Pairs after NaNs and 1s removal,', '# Pairs after area cut,', 'Partial Corr,','Rev. Partial Corr','\n']);
fprintf(fpCSV3,['RunNr,','Run,','Group,','# Cells,','# Pairs after cuts,', ...
                'RegSlope S1-B,', '95ConfIntLB S1-B,', '95ConfIntUB S1-B,','R2 S1-B,','F-value S1-B,','p-value S1-B,', ...
                'RegSlope S2-B,', '95ConfIntLB S2-B,', '95ConfIntUB S2-B,','R2 S2-B,','F-value S2-B,','p-value S2-B,', ...
                '\n']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Set RUN PARAMTERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% run parameters
runQuiet = 1                % run verbose (0) or quiet (1)

%set limits for cuts on cell firing properties
SpikeCountLowerLimit_S1 = 30 % sleep 1 lower spike count limit
SpikeCountLowerLimit_S2 = 30 % sleep 2 lower spike count limit
SpikeCountLowerLimit_B  = 90 % behavior epoch lower spike count limit
rateLim = Inf                % mean firing rate upper limit (Hz) for cells during behavior for inclusion in the analysis

%set restriction intervals for sleep 1 and 2 epochs
sleep1_offset_sec = 60      % 60 sec offset before end of sleep 1. 
sleep1_interval_sec = -900  % restrict sleep1 to 10 min interval before end of sleep1 minus offset
sleep2_offset_sec = 180     % give animal 3 min to start sleep2 
sleep2_interval_sec = 900   % restrict sleep2 to 10 min interval after offset

% paramters for CCG-matrix and generation
exclude_sameTT = 1          % exclude cell pairs from same tetrode  
binsize_S = 10              % sleep binsize:  5 for 100ms windows or 2.5 for 50 ms windows
binsize_B = 10              % behavior binsize: binsize for CCGs in ms
numbins = 100               % number of bins (also timelags) for CCGs
normal = 1                  % if normal = 1, correct for outliers in bias (often resulting in 1 or -1 for bias) 

% cut on area of crosscorrelogram of cell pairs
area_cut = 1
%area_threshold = 0         % Run xx0
%area_threshold = exp(8)    % Run xx1
%area_threshold = exp(10)    % Run xx2
area_threshold = [0, exp(7:0.5:12)]    % Run 010 030

% load the session list
sessList = SessionList02;

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

NRuns = length(area_threshold);
for irun = 1:NRuns
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
	S1_bias_pooled = [];
	B_bias_pooled = [];
	S2_bias_pooled = [];
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
        
        clear tfiles_*;
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
            if ~runQuiet
                fprintf('\n RAT %s, Session %s:\n session has N = %d cells prior to cuts\n\n',ses.animal,ses.name,NCellsPreCuts); 
            end            
            
            SS1 = getSpikeMatrix(ses,'sleep1');
            SB  = getSpikeMatrix(ses,'maze1');
            SS2 = getSpikeMatrix(ses,'sleep2');
            
            sameTT = getSameNTrodeMatrix(ses);
            
            % restrict the S-matrices to the desired sub-intervals:
            SS1 = RestrictTheSMatrix2Interval(SS1,sleep1_offset_sec,sleep1_interval_sec);
            SS2 = RestrictTheSMatrix2Interval(SS2,sleep2_offset_sec,sleep2_interval_sec);
            
            % check intervals
            [ts0_S1,ts1_S1] = SessionStartEndTS(SS1);
            SS1_interval_sec = (ts1_S1-ts0_S1)/10000;
            [ts0_B,ts1_B] = SessionStartEndTS(SB);
            SB_interval_sec = (ts1_B-ts0_B)/10000;
            [ts0_S2,ts1_S2] = SessionStartEndTS(SS2);
            SS2_interval_sec = (ts1_S2-ts0_S2)/10000;
            if ~runQuiet            
                fprintf('SS1: ts0 = %d, ts1 = %d, interval(sec) = %f\n',ts0_S1,ts1_S1,SS1_interval_sec);
                fprintf('SB: ts0 = %d, ts1 = %d, interval(sec) = %f\n',ts0_B,ts1_B,SB_interval_sec);
                fprintf('SS2: ts0 = %d, ts1 = %d, interval(sec) = %f\n\n',ts0_S2,ts1_S2,SS2_interval_sec);
            end
            
            % check avg firing rates
            avg_rate_S1 = SMatrix_getAvgFiringRates(SS1);
            avg_rate_B = SMatrix_getAvgFiringRates(SB);
            avg_rate_S2 = SMatrix_getAvgFiringRates(SS2);
            if ~runQuiet
                fprintf('Average Firing Rates prior to cuts:\n');
                fprintf('SS1: avgRate = %f +- std = %f (Hz)\n',mean(avg_rate_S1),std(avg_rate_S1));
                fprintf('SB:  avgRate = %f +- std = %f (Hz)\n',mean(avg_rate_B),std(avg_rate_B));
                fprintf('SS2: avgRate = %f +- std = %f (Hz)\n\n',mean(avg_rate_S2),std(avg_rate_S2));
            end
            
            % apply cuts on tfiles: compile a list of cells to be discarded for this analysis
            if ~runQuiet
                disp('applying cuts on tfiles:');
            end
            badindx = [];        
            for i = 1:length(SS1)
                secs = (EndTime(SB{i})-StartTime(SB{i}))/10000;
                
                if (length(data(SS1{i})) < SpikeCountLowerLimit_S1)     % must be at least 30 spikes during S1
                    badindx = [ badindx ; i ];
                    if ~runQuiet
                        fprintf('The tfile %s had less than %d spikes during S1 \n',tfiles_S1{isess}{i},SpikeCountLowerLimit_S1);
                    end
                    
                elseif (length(data(SS2{i})) < SpikeCountLowerLimit_S2)     % must be at least 30 spikes during S2
                    badindx = [ badindx ; i ];
                    if ~runQuiet
                        fprintf('The tfile %s had less than %d spikes during S2 \n',tfiles_S2{isess}{i},SpikeCountLowerLimit_S2);
                    end
                    
                elseif (length(data(SB{i})) < SpikeCountLowerLimit_B)        %must be at least 90 spikes during behavior
                    badindx = [ badindx ; i];
                    if ~runQuiet
                        fprintf('The tfile %s had less than %d spikes during Behavior \n',tfiles_B{isess}{i},SpikeCountLowerLimit_B);
                    end
                    
                elseif ((length(data(SB{i}))/secs) > rateLim)  % Check mean firing rate during behavior
                    badindx = [ badindx ; i];
                    if ~runQuiet
                        fprintf('The tfile %s had a mean rate greater than 2.2 Hz during behavior \n',tfiles_B{isess}{i});
                    end        
                end %if
            end
            
            % discard the list of cells not passing above criteria
            if ~isempty(badindx)
                if ~runQuiet
                    disp('Found tfiles that did not pass all cuts')
                    fprintf('%d cells will be removed from the analysis.\n',length(badindx));   
                end
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
            if ~runQuiet
                fprintf('\n RAT %s, Session %s:\n session has N = %d cells and nPairs = %d valid pairs after cuts\n\n',ses.animal,ses.name,length(SB),nPairs); 
            end
            
            if nPairs == 0
                % ignore current session
                if ~runQuiet
                    fprintf('\n ===> RAT %s, Session %s:\n this session has not enough pairs (Pairs = %d) and is skipped!!\n\n',ses.animal,ses.name,nPairs); 
                end
                continue;
            end
            
            % check avg firing rates again after cuts
            avg_rate_S1 = SMatrix_getAvgFiringRates(SS1);
            avgRateSS1{irat,isess} = [mean(avg_rate_S1),std(avg_rate_S1)];
            avg_rate_B = SMatrix_getAvgFiringRates(SB);
            avgRateSB{irat,isess} = [mean(avg_rate_B),std(avg_rate_B)];
            avg_rate_S2 = SMatrix_getAvgFiringRates(SS2);
            avgRateSS2{irat,isess} = [mean(avg_rate_S2),std(avg_rate_S2)];
            if ~runQuiet
                fprintf('Average Firing Rates after cuts:\n');
                fprintf('SS1: avgRate = %f +- std = %f (Hz)\n',mean(avg_rate_S1),std(avg_rate_S1));
                fprintf('SB:  avgRate = %f +- std = %f (Hz)\n',mean(avg_rate_B),std(avg_rate_B));
                fprintf('SS2: avgRate = %f +- std = %f (Hz)\n\n',mean(avg_rate_S2),std(avg_rate_S2));
            end
            
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
            fprintf(fpCSV1,['%d,','%s,'     ,'%s,'   ,'%s,','%f,'           ,'%f,'          ,'%f,'           ,'%d,'        ,'%d,','%d'           ,'\n'],...
                            irun,ses.animal,ses.name,group,SS1_interval_sec,SB_interval_sec,SS2_interval_sec,NCellsPreCuts, N   , NParisPreCuts );
            
        end % for isess
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End loop over sessions
        
        if nsess == 0
            % ignore current animal
            fprintf('\n ===> run %d, RAT %s has contributed 0 sessions to pool and is skipped!!\n\n',irun,ses.animal); 
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
        fprintf('\n Run %d: Removing %d pairs (rows) with NaNs from the Bias Matrix, leaving %d pairs.' ,irun,length(ix_rows),nr)
        if (normal)              
            % Remove all of the -1 or 1 values from the Bias, this usually occurs
            % when only one or two bins on one side has any values such that, e.g.
            % the sum of -200ms to 0ms <<<<< sum of 0ms to 200ms.
            ix_rows = find(sum(abs(ABM(:,1:3))==1,2));
            ABM(ix_rows,:) = [];
            preM(ix_rows,:) = [];
            postM(ix_rows,:) = [];
            [nr, nc] = size(ABM);
            fprintf('\n Run %d: Removing %d pairs (rows) with 1s from the Bias Matrix, leaving %d pairs.' ,irun,length(ix_rows),nr)
        end %if             
        NPairsPostNaNsCut = nr;
        
        % AREA CUT:
        % Remove all rows from the Bias matrix which correspond to pairs with area in Crosscorrlogram < area_cut
        if area_cut
            areaM = preM+postM;
            ix_rows = find(sum(areaM < area_threshold(irun),2));
            areaM(ix_rows,:) = [];
            ABM(ix_rows,:) = [];
            preM(ix_rows,:) = [];
            postM(ix_rows,:) = [];
            [nr, nc] = size(ABM);
            fprintf('\n Run %d: Removing %d pairs (rows) below area_threshold %f from the Bias Matrix, leaving %d pairs.' ,irun,length(ix_rows),area_threshold(irun),nr)
        end %if
        NPairsPostAreaCut = nr;
        
        % accumulate areaM for all rats
        areaMall = [areaMall; areaM];
        
        
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
        fprintf('\n Run %d, RAT %s, %d sessions pooled together into %d cells:\n EV = %f, EV_rev = %f \n',...
            irun, ses.animal, nsess, NCellsPerRat, S2B_allbiaspart(1), S2B_allbiaspart_rev(1)); 
        fprintf('\n=============================================================================================\n');
                
        animalSessionList = sessListPerRat{irat};
        NCellsTotal = NCellsTotal + NCellsPerRat;
        
        % save rat info to CSV
        fprintf(fpCSV2,['%d,','%s,','%s,','%d,','%d,'        ,'%d,'         , '%d,'            , '%d,'            , '%f,'             ,'%f'                ,'\n'], ...
                   irun,ses.animal,group,nsess, NCellsPerRat, NPairsPreCuts, NPairsPostNaNsCut, NPairsPostAreaCut, S2B_allbiaspart(1), S2B_allbiaspart_rev(1));
        
        % save the data to the a file that is called "RAT#_bias_results.mat"
        %filename = [outputDir, filesep, 'Run',runNrStr,'_', 'Bias_results_', ratlist{irat}, '_', num2str(binsize_S *20), 'ms'];
        %fprintf('\n writing variables to file %s in dir %s' ,filename,pwd)
        %save(filename,'bias*','pre*','post*','CCG*','*allbia*','*sessbia*','ABM','NCellsPerRat','animalSessionList','nsess');
        
        % accumulate ABMs for regression 
        S1_bias_pooled = [S1_bias_pooled; ABM(:,1)];
        B_bias_pooled = [B_bias_pooled; ABM(:,2)];
        S2_bias_pooled = [S2_bias_pooled; ABM(:,3)];
        
        
        % Clear out the memory and move on to the next rat
%         clear *sessbia*; 
%         clear *allbia*; 
%         clear CCG*; 
%         clear S1* S2* SS* SB*; 
%         clear sameTT;
%         clear sys*; 
%         clear tfiles*; 
%         clear ts*; 
        %pack;
        
        
    end% for irat
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End loop over rats
    
    
    % Regression Slopes
    [NDots,Ncols] = size(S1_bias_pooled);
    X = [ones(size(B_bias_pooled,1),1) B_bias_pooled];
    if NDots > 2
		[S1Bstats.b,S1Bstats.bint,S1Bstats.r,S1Bstats.rint,S1Bstats.stats] = regress(S1_bias_pooled,X);  % S1-B 
		
		fh = figure;
		scrPos = get(0,'ScreenSize');
		set(fh,'Units','pixels');
		set(fh,'Position',[0,0,(scrPos(4)-20)/2,scrPos(4)-20]);
		orient tall;
		subplot(2,1,1);
		plot(B_bias_pooled,S1_bias_pooled,'.') 
		hold on
		reg = S1Bstats;
		x = [-1.5,1.5];
		y = reg.b(2)*x + reg.b(1);
		lh = line(x,y);
		set(lh,'LineStyle','-','Color','red');
		textstr = sprintf('%6.3f*x %+6.3f',reg.b(2),reg.b(1));
		text(-1,0.9,textstr);
		hold off
		grid on
		title(['Run' runNrStr ' ' num2str(irun) ' ' group ': ' num2str(length(ratlist)) 'Rats combined: B vs S1    (' num2str(NDots) ' data points)'])
		xlabel('Behavior')
		ylabel('Sleep 1')
		axis([-1,1,-1,1])
		axis equal
	
	
		[S2Bstats.b,S2Bstats.bint,S2Bstats.r,S2Bstats.rint,S2Bstats.stats] = regress(S2_bias_pooled,X);  % S2-B 
		
		subplot(2,1,2);
		plot(B_bias_pooled,S2_bias_pooled,'.') 
		hold on
		reg = S2Bstats;
		x = [-1.5,1.5];
		y = reg.b(2)*x + reg.b(1);
		lh = line(x,y);
		set(lh,'LineStyle','-','Color','red');
		textstr = sprintf('%6.3f*x %+6.3f',reg.b(2),reg.b(1));
		text(-1,0.9,textstr);
		hold off
		grid on
		title(['Run' runNrStr ' ' num2str(irun) ' ' group ': ' num2str(length(ratlist)) ' Rats combined: B vs S2    (' num2str(NDots) ' data points)'])
		xlabel('Behavior')
		ylabel('Sleep 2')
		axis([-1,1,-1,1])
		axis equal
		
		saveas(gcf,[outputDir, filesep, 'Run',runNrStr,'_',num2str(irun),'_',  outPrefix, group,'_AllRatsCombined'],'fig');

    else
        S1Bstats.b = [NaN,NaN];
        S1Bstats.bint = [NaN,NaN;NaN,NaN];
        S1Bstats.stats = [NaN,NaN,NaN];
        S2Bstats.b = [NaN,NaN];
        S2Bstats.bint = [NaN,NaN;NaN,NaN];
        S2Bstats.stats = [NaN,NaN,NaN];
    end
    
    
    % write regression result to CVS
    fprintf(fpCSV1,'\n')
    fprintf(fpCSV2,'\n')    
    %['RunNr','Run','Group','# Cells','# Pairs after cuts',...
    %       'Regression Slope S1-B', '95% ConfInterval Lower Bound  S1-B', '95% ConfInterval Upper Bound S1-B','R^2 S1-B','F-value S1-B','p-value S1-B', ...
    %       'Regression Slope S2-B', '95% ConfInterval Lower Bound  S2-B', '95% ConfInterval Upper Bound S2-B','R^2 S2-B','F-value S2-B','p-value S2-B', ...
    %       '\n']
    fprintf(fpCSV3,...
        ['%d,','%d,','%s,','%d,'      ,'%d,' , ...
        '%f,'        , '%f,'             , '%f,'             ,'%f,'             ,'%f,'            ,'%f,'         , ...
        '%f,'        , '%f,'             , '%f,'             ,'%f,'             ,'%f,'            ,'%f,'         , '\n'], ...
        runNr,irun,group,NCellsTotal, NDots, ...
        S1Bstats.b(2), S1Bstats.bint(2,1), S1Bstats.bint(2,2), S1Bstats.stats(1),S1Bstats.stats(2),S1Bstats.stats(3), ... 
        S2Bstats.b(2), S2Bstats.bint(2,1), S2Bstats.bint(2,2), S2Bstats.stats(1),S2Bstats.stats(2),S2Bstats.stats(3));
    
end% for irun
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End of loop over runs with different parameters

% save figure
%filename = [outputDir filesep outPrefix 'BiasHistograms_AllAnimals_AllSessions' num2str(Btime_binsize) '_' limname 'Hz_Results'];
%saveas(fh,filename,'fig');

fclose('all');
diary off

