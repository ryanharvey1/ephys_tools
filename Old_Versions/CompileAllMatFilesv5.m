% CompileAllMatFilesv5
% This script will cycle through a parent directory and pull out info from the .mat files
% input:
%       Laura=0; If you want Laura's data
%       Ryan=1; If you want Ryan's data
% output:
%       SAVES EXCEL SHEETS FOR EACH RAT WITH ALL VARIABLES
%
% Ryan Harvey 9/16
function CompileAllMatFilesv5(Laura,Ryan)
tic
addpath(genpath('F:\Users\BClarkLab\GoogleDrive\MatlabDir\BClarkToolbox\Analyses\spikeCode'));
addpath(genpath('d:\Users\BClarkLab\GoogleDrive\MatlabDir\BClarkToolbox\Analyses\spikeCode'));

if Ryan==1
    rats={'RH11','RH13','RH14','RH16','LS17','LS19','LS21','LS23','LE2813','LE2821'};
    datafolder='F:\Users\reharvey\Place_Cell_Data\PAE_Rat\';
elseif Laura==1
    rats={'LB01','LB02','LB03','LB04','LB05','LB06'};
    datafolder='D:\ClarkP30_Recordings\';
end

for irats=1:length(rats)
    parent=strcat(datafolder,rats(irats));
    disp(['CYCLING THROUGH RAT:',char(rats(irats))])
    parent=char(parent);
    structdir=dir(parent);
    K=1;
    for I=1:length(structdir) % 1 TO # OF FILES IN DIR
        if structdir(I).isdir && structdir(I).name(1) ~= '.' % IF DIR NOT '.'
            if exist(([parent filesep structdir(I).name filesep 'TT']),'file'); % LOCATE TT FOLDER
                cd([parent filesep structdir(I).name filesep 'TT']); % CD TO TT FOLDER
                CurentFolder=pwd;
                CurrentMat=dir('*.mat'); % LOCATE .MAT FILES IN TT FOLDER
                for J=1:length(CurrentMat) % 1 TO # OF .MAT FILES IN TT FOLDER
                    CurrentMatworking=CurrentMat(J).name;
                    if sum(ismember(CurrentMatworking,'spikeData'))==11 && sum(ismember(CurrentMatworking,'pathProperties'))~=16
                        disp(['READING: ',CurentFolder filesep CurrentMatworking])
                        warning('off','MATLAB:load:variableNotFound');
                        load(CurrentMatworking,'linear_track'); if isequal(linear_track,'no'); continue; end
                        load(CurrentMatworking,'linear_track','InformationContent','sparsity','Sparsity','PeakRate','OverallFR','DistanceFromTrackEnd',...
                            'FieldWidth','mean_vector_length','preferred_Direction','Direct_infoContent','Stats','nSpikes','DirectionalityIndex',...
                            'Coherence','BasicLoco','ThPrecess','ThetaStats','PerSpkslessthan2ms','Displacement','infieldFR','outfieldFR','CI'); % OPEN .MAT WITH THESE VARS
                        if ~exist('Displacement','var') || ~exist('Stats','var');
                            continue; 
                        end
                        if ~exist('ThPrecess','var'); % ONLY LINEAR TRACK HAS PH PRECESSION, CURRENTLY
                            ThPrecess.Correlation=NaN;ThPrecess.DOM=NaN;ThPrecess.slope=NaN;ThPrecess.meanFR=NaN; ThPrecess.RSquared=NaN;
                        end
                        if ~exist('PerSpkslessthan2ms','var'); PerSpkslessthan2ms=NaN; end
                        if ~exist('ThetaStats','var')
                            ThetaStats.MeanThetaFreq=NaN;ThetaStats.MeanOverallPow=NaN;ThetaStats.ThetaRatio=NaN;
                        end
                        if ~exist('DirectionalityIndex','var'); DirectionalityIndex=[];end;
                        if ~exist('BasicLoco','var')==1; BasicLoco.AverageAnglePerSec=NaN; BasicLoco.OverallDistTraveled=NaN; BasicLoco.MeanVelocity=NaN; end
%                         if exist('sparsity','var')==1 
%                             [TTclust]= strsplit(CurrentMatworking,'_');
%                             if exist(char(strcat(TTclust(1), '_ClusterSummary_', TTclust(2),'.mat')),'file')>0
%                                 load(char(strcat(TTclust(1), '_ClusterSummary_', TTclust(2),'.mat')),'CI');
                            if exist('CI','var')==0
                                CI.CluSep.IsolationDist=NaN; CI.CluSep.L_Ratio=NaN;
                                CI.PeakToTroughPts=[NaN,NaN,NaN,NaN]; CI.PeakHalfWidthPts=[NaN,NaN,NaN,NaN];
                            end
                                if isfield(CI.CluSep,'IsolationDist')==0; CI.CluSep.IsolationDist=NaN; CI.CluSep.L_Ratio=NaN; end
                                if isfield(CI,'PeakToTroughPts')==0; CI.PeakToTroughPts=[NaN,NaN,NaN,NaN]; CI.PeakHalfWidthPts=[NaN,NaN,NaN,NaN]; end

%                         elseif exist('Sparsity','var')==1
%                             continue
%                         end
                        if exist('infieldFR','var')==1 % IF INFO CONTENT EXISTS IN A GIVEN .MAT FILE, IT'S PROBABLY THE RIGHT FILE
                            % COMBINE INTO CELL ARRAY
                            try
                            AllSpikeData(K,:)=[[CurentFolder filesep CurrentMatworking],num2cell(InformationContent),num2cell(Coherence),num2cell(sparsity),num2cell(PeakRate),...
                                num2cell(OverallFR),num2cell(max(DistanceFromTrackEnd)),num2cell(FieldWidth),...
                                num2cell(CI.PeakToTroughPts(1)),num2cell(CI.PeakToTroughPts(2)),num2cell(CI.PeakToTroughPts(3)),num2cell(CI.PeakToTroughPts(4)),...
                                num2cell(CI.PeakHalfWidthPts(1)),num2cell(CI.PeakHalfWidthPts(2)),num2cell(CI.PeakHalfWidthPts(3)),num2cell(CI.PeakHalfWidthPts(4)),...
                                num2cell(nSpikes), num2cell(CI.CluSep.IsolationDist), num2cell(CI.CluSep.L_Ratio),num2cell(mean_vector_length),...
                                num2cell(preferred_Direction(1)),num2cell(Direct_infoContent),num2cell(DirectionalityIndex),...
                                num2cell(Stats.MeanResultantPhase(1)),num2cell(Stats.MeanResultantPhase(2)),num2cell(Stats.MeanResultantPhase(3)),...
                                num2cell(Stats.MeanResultantPhase(4)),num2cell(Stats.MeanResultantPhase(5)),num2cell(Stats.MeanResultantPhase(6)),...
                                num2cell(Stats.RayleighsTest.PVal(1)),num2cell(Stats.RayleighsTest.PVal(2)),num2cell(Stats.RayleighsTest.PVal(3)),...
                                num2cell(Stats.RayleighsTest.PVal(4)),num2cell(Stats.RayleighsTest.PVal(5)),num2cell(Stats.RayleighsTest.PVal(6)),...
                                num2cell(Stats.RayleighsTest.ZValue(1)),num2cell(Stats.RayleighsTest.ZValue(2)),num2cell(Stats.RayleighsTest.ZValue(3)),...
                                num2cell(Stats.RayleighsTest.ZValue(4)),num2cell(Stats.RayleighsTest.ZValue(5)),num2cell(Stats.RayleighsTest.ZValue(6)),...
                                num2cell(Stats.Meanz(1)),num2cell(Stats.Meanz(2)),num2cell(Stats.Meanz(3)),...
                                num2cell(Stats.Meanz(4)),num2cell(Stats.Meanz(5)),num2cell(Stats.Meanz(6)),...
                                num2cell(Stats.skewness(1)),num2cell(Stats.skewness(2)),num2cell(Stats.skewness(3)),...
                                num2cell(Stats.skewness(4)),num2cell(Stats.skewness(5)),num2cell(Stats.skewness(6)),...
                                num2cell(Stats.kurtosis(1)),num2cell(Stats.kurtosis(2)),num2cell(Stats.kurtosis(3)),...
                                num2cell(Stats.kurtosis(4)),num2cell(Stats.kurtosis(5)),num2cell(Stats.kurtosis(6)),...
                                num2cell(Stats.median(1)),num2cell(Stats.median(2)),num2cell(Stats.median(3)),...
                                num2cell(Stats.median(4)),num2cell(Stats.median(5)),num2cell(Stats.median(6)),...
                                num2cell(BasicLoco.AverageAnglePerSec),num2cell(BasicLoco.OverallDistTraveled),num2cell(BasicLoco.MeanVelocity),...
                                num2cell(ThPrecess.Correlation),num2cell(ThPrecess.DOM),num2cell(ThPrecess.slope),num2cell(ThPrecess.meanFR),num2cell(ThPrecess.RSquared),...
                                num2cell(ThetaStats.MeanThetaFreq),num2cell(ThetaStats.MeanOverallPow),num2cell(ThetaStats.ThetaRatio),num2cell(PerSpkslessthan2ms),num2cell(Displacement),...
                                num2cell(infieldFR),num2cell(outfieldFR)]; 
                            catch
                                test=1
                            end
                            K=K+1; % CYCLE THROUGH OUTPUT FILE SO YOU DON'T OVERWRITE
                        end
                        if ~exist('AllSpikeData','var') % IF THE .MAT FILE WAS INCORRECT
                            keep('I','J','K','structdir','parent','datafolder','CurrentMat','CurrentFolder','rats','irats');
                        else
                            keep('I','J','K','AllSpikeData','structdir','parent','datafolder','CurrentMat','CurentFolder','rats','irats'); % KEEP ALLSPIKEDATA & STRUCTDIR IF .MAT WAS CORRECT
                        end
                    end
                end
                if exist('AllSpikeData','var')
                    keep('I','K','AllSpikeData','structdir','parent','datafolder','CurrentMat','rats','irats'); % DON'T KEEP VAR J BECAUSE YOU JUST EXITED THE J LOOP
                else
                    keep('I','K','structdir','parent','datafolder','CurrentMat','rats','irats'); 
                end
            end
        end
    end
    if exist('AllSpikeData','var')
        AllSpikeData = cell2table(AllSpikeData,'VariableNames',{'CurentDir','InformationContent','Coherence','Sparsity','PeakRate','OverallFiringRate','DistanceFromTrackEnd','FieldWidth',...
            'PeaktoTrough1','PeaktoTrough2','PeaktoTrough3','PeaktoTrough4','PeakHalfWidthPts1','PeakHalfWidthPts2','PeakHalfWidthPts3','PeakHalfWidthPts4',...
            'nSpikes','IsolationDist','L_Ratio','mean_vector_length','preferred_Direction','Direct_infoContent','DirectionalityIndex',...
            'RLength_Delta','RLength_Theta','RLength_Alpha','RLength_Beta','RLength_Gamma','RLength_HighGamma',...
            'Rayleigh_Delta','Rayleigh_Theta','Rayleigh_Alpha','Rayleigh_Beta','Rayleigh_Gamma','Rayleigh_HighGamma',...
            'RayleighZ_Delta','RayleighZ_Theta','RayleighZ_Alpha','RayleighZ_Beta','RayleighZ_Gamma','RayleighZ_HighGamma',...
            'MeanPhaseDeg_Delta','MeanPhaseDeg_Theta','MeanPhaseDeg_Alpha','MeanPhaseDeg_Beta','MeanPhaseDeg_Gamma','MeanPhaseDeg_HighGamma',...
            'skewness_Delta','skewness_Theta','skewness_Alpha','skewness_Beta','skewness_Gamma','skewness_HighGamma',...
            'kurtosis_Delta','kurtosis_Theta','kurtosis_Alpha','kurtosis_Beta','kurtosis_Gamma','kurtosis_HighGamma',...
            'median_Delta','median_Theta','median_Alpha','median_Beta','median_Gamma','median_HighGamma'...
            'AverageAnglePerSec','OverallDistTraveled','MeanVelocity',...
            'PhPrecessCorr','DepthofModulation','PhPrecessSlope','PhPrecessMeanFR','PH_RSquared',...
            'MeanThetaFreq','MeanOverallPow','ThetaRatio','PerSpkslessthan2ms','Displacement','infieldFR','outfieldFR'}); % CREATE TABLE
        [filepath, filename] = fileparts(parent);
        if ismac==1
            writetable(AllSpikeData,[filepath filesep '_AllSpikeData' filename '.csv'],'WriteRowNames',true);
        else
            delete([filepath filesep '_AllSpikeData' filename '.xlsx'])
            writetable(AllSpikeData,[filepath filesep '_AllSpikeData' filename '.xlsx']);
        end
    end
    disp(['DONE WITH:',char(rats(irats))])
    keep('rats','irats','parent','datafolder')
end
disp('done')
toc




