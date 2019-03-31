%% deleteJPEGS
% This script will cycle through a parent directory and pull out info from the .mat files
% Ryan Harvey 12/16
clear,close all
if ismac==0
    addpath(genpath('F:\Users\BClarkLab\GoogleDrive\MatlabDir\BClarkToolbox\Analyses\spikeCode'));
    addpath(genpath('d:\Users\BClarkLab\GoogleDrive\MatlabDir\BClarkToolbox\Analyses\spikeCode'));
    addpath(genpath('F:\Users\BClarkLab\GoogleDrive\MatlabDir\BClarkToolbox\Analysis'));
    
    for igroup=1:2
        SmoothRateMap_Rightw=[];
        SmoothRateMap_Leftw=[];
%         storingFile=[];
        if igroup==1
            rats={'RH13','RH14','LS21','LS23','LE2821'};
        elseif igroup==2
            rats={'RH11','RH16','LS17','LS19','LE2813'};
        end
        for irats=1:length(rats)
            parent=strcat('F:\Users\reharvey\Place_Cell_Data\PAE_Rat\',rats(irats));
            disp(['CYCLING THROUGH RAT:',char(rats(irats))])
            parent=char(parent);
            structdir=dir(parent);
                    K=1;
            for I=1:length(structdir) % 1 TO # OF FILES IN DIR
                if structdir(I).isdir && structdir(I).name(1) ~= '.' % IF DIR NOT '.'
                    if exist(([parent filesep structdir(I).name filesep 'TT']),'file'); % LOCATE TT FOLDER
                        CurentFolder=cd([parent filesep structdir(I).name filesep 'TT']); % CD TO TT FOLDER
                        CurrentMat=dir('*.mat'); % LOCATE .MAT FILES IN TT FOLDER
                        for J=1:length(CurrentMat) % 1 TO # OF .MAT FILES IN TT FOLDER
                            CurrentMatworking=CurrentMat(J).name;
                            if sum(ismember(CurrentMatworking,'spikeData'))==11 && sum(ismember(CurrentMatworking,'pathProperties'))~=16
                                warning('off','MATLAB:load:variableNotFound');
                                load(CurrentMatworking,'linear_track'); if isequal(linear_track,'no'); continue; end
                                load(CurrentMatworking,'right','left','iruns','InformationContent','PeakRate','OverallFR','nSpikes'); % OPEN .MAT WITH THESE VARS
                                disp(['READING: ',CurentFolder filesep CurrentMatworking])
                                if exist('InformationContent','var') && exist('right','var')
                                    if InformationContent>=0.30 && nSpikes>=50 && PeakRate>=.25 %&& OverallFR<=10
                                        % arrange ratemaps so that SmoothRateMap_Right aways contains the place field
                                        if iruns==1
                                            SmoothRateMap_Right=right.SmoothRateMap; SmoothRateMap_Left=left.SmoothRateMap;
                                        elseif iruns==2
                                            SmoothRateMap_Right=left.SmoothRateMap; SmoothRateMap_Left=right.SmoothRateMap;
                                        end
                                        % resize to 22 bins and concat
                                        if size(SmoothRateMap_Right,2)<37; SmoothRateMap_Right=imresize(SmoothRateMap_Right,[1,37]);end
                                        SmoothRateMap_Rightw=[SmoothRateMap_Rightw;SmoothRateMap_Right];
                                        
                                        if size(SmoothRateMap_Left,2)<37; SmoothRateMap_Left=imresize(SmoothRateMap_Left,[1,37]);end
                                        SmoothRateMap_Leftw=[SmoothRateMap_Leftw;SmoothRateMap_Left];
                                        storingFile{K,:}=[CurentFolder filesep CurrentMatworking]; K=K+1;
                                    end
                                end
                                if ~exist('SmoothRateMap_Right','var') % IF THE .MAT FILE WAS INCORRECT
                                    keep('storingFile','I','J','K','structdir','parent','CurrentMat','CurrentFolder','rats','irats','igroup');
                                else
                                    keep('storingFile','I','J','K','SmoothRateMap_Rightw','SmoothRateMap_Leftw','structdir','parent','CurrentMat','CurentFolder','rats','irats','igroup'); % KEEP ALLSPIKEDATA & STRUCTDIR IF .MAT WAS CORRECT
                                end
                            end
                        end
                        keep('storingFile','I','K','SmoothRateMap_Rightw','SmoothRateMap_Leftw','structdir','parent','CurrentMat','rats','irats','igroup'); % DON'T KEEP VAR J BECAUSE YOU JUST EXITED THE J LOOP
                    end
                end
            end
            disp(['DONE WITH:',char(rats(irats))])
            keep('storingFile','rats','irats','parent','igroup','SmoothRateMap_Rightw','SmoothRateMap_Leftw')
        end
        
        
        % WORKING******** USE storingFile TO DELETE MULTI-PAIRS
        
        
%         SmoothRateMap_Rightw=unique(SmoothRateMap_Rightw,'rows','stable');
%         SmoothRateMap_Leftw=unique(SmoothRateMap_Leftw,'rows','stable');
        
        if igroup==1
            Group=('Control');
        elseif igroup==2
            Group=('PAE');
        end
        % TAKE FILTERED COMPILED DATA AND PLOT THEM FOR EACH ANIMAL IN EACH DIRECTION
        DirectionalMatrix(parent,[],0,SmoothRateMap_Rightw,SmoothRateMap_Leftw,Group);
        % CREATE POPULATION MAPS
        %     DirectionalMatrix(parent,[],0,PopSmoothRateMap,Group);
        keep('igroup','parent')
    end
    % MAKE FIGURES AND COMPUTE STATS
    %
    [filepath, ~] = fileparts(parent);
else
    addpath(genpath('/Users/ryanharvey/GoogleDrive/MatlabDir/NSMA_Toolbox'))
    addpath(genpath('F:\Users\BClarkLab\GoogleDrive\MatlabDir\NSMA_Toolbox'))
    for igroup=1:2
        if igroup==1
            Group=('Control');
        elseif igroup==2
            Group=('PAE');
        end
        if ismac==1
            load(['/Users/RyanHarvey/Dropbox/school work/UNM/Lab/Projects/PAE Project/PLACECELL_DATASHEETS' filesep Group '_Data.mat'])
        else
            load([filepath filesep Group '_Data.mat'])
        end
        
        % REMOVE NANS
        SmoothRateMap_Right_arrangedR(isnan(SmoothRateMap_Right_arrangedR))=0;
        SmoothRateMap_Left_arrangedR(isnan(SmoothRateMap_Left_arrangedR))=0;
        
        % AUTO CORRELATE & CROSS CORRELATE
        autoR=corrcoef_AB(SmoothRateMap_Right_arrangedR,SmoothRateMap_Right_arrangedR);
        R=corrcoef_AB(SmoothRateMap_Right_arrangedR,SmoothRateMap_Left_arrangedR);
        
%         [M,I]=max(R);
        
        % Find sig correlations
        [~,P]=corr(SmoothRateMap_Right_arrangedR,SmoothRateMap_Left_arrangedR);
        P(P>.0001)=NaN;
        maxP=max(max(P));
        P(isnan(P))=maxP;
        
        % ROW BY ROW CORRELATION 
        RowCorrelations=zeros(size(SmoothRateMap_Right_arrangedR,1),1);
        for irow=1:size(SmoothRateMap_Right_arrangedR,1)
            RowCorrelations(irow,1)=corr2(SmoothRateMap_Right_arrangedR(irow,:),SmoothRateMap_Left_arrangedR(irow,:));
        end
         
        [f,x] = ecdf(RowCorrelations);
        if igroup==1 
            cdfControl=[f,x]; RowCorrelations_control=RowCorrelations;
        else
            cdfPAE=[f,x]; RowCorrelations_PAE=RowCorrelations;
        end

        % FIND DECORRELATING DISTANCE
        decordist=zeros(size(R,2),1);
        decordist2=zeros(size(R,2),1);
        for a=0:size(R,2)-1
            decordist(a+1)=nanmean(diag(R, a));
            decordist2(a+1)=nanmean(diag(R, -a));
        end
        reversedecorr=flipud(decordist2); reversedecorr(end)=[];
        decordistALL(:,igroup)=[reversedecorr;decordist];
        
        if igroup==1 
            diagCorr_control=diag(R, 0);
        else
            diagCorr_PAE=diag(R, 0);
        end
        
        %PLOT POPULATION VECTORS
        figure(5); subplot(3,2,igroup);
        matfordis=[autoR,R;R,autoR];
        h = pcolor(matfordis);
        colormap jet
        axis off
        colorbar
        box off
        set(h, 'EdgeColor', 'none');
        set(gca,'YDir','reverse');
        title(['Population vector correlation: ',num2str(Group)]);
        
        %SIG PVALUES OF POPULATION VECTOR
        subplot(3,2,igroup+2); 
        p2=pcolor(-P);
        colormap jet
        axis off
        box off
        set(p2,'EdgeColor','none');
        set(gca,'YDir','reverse');
        title(['Significant Population Correlations : ',num2str(Group)]);
    end
    % CUMULATIVE FREQ OF ROW BY ROW CORRELATIONS
    subplot(3,2,[5.4 5.6])
    p3=plot(cdfControl(:,2),cdfControl(:,1),'DisplayName','Control');
    set(p3,'LineWidth',4,'Color','k')
    hold on
    p3_2=plot(cdfPAE(:,2),cdfPAE(:,1));
    set(p3_2,'LineWidth',4,'Color','r','DisplayName','PAE')
    box off
    legend('show','Location','best')
    xlabel('Directional Correlation');ylabel('Cumulative Frequency')
    
    %         [p,h,stats] = signrank(cdfControl(:,2),cdfPAE(:,2))
    
    % Mean Half-width
    half1=decordistALL(:,1);
    half2=decordistALL(:,2);
    halfwidth1=median(half1(half1>0));
    halfwidth2=median(half2(half2>0));
    
    x=linspace(-120,120,73);
    x1=interp1(half1(1:37,1),x(1,1:37),halfwidth1,'linear');
    x2=interp1(half1(37:end,1),x(1,37:end),halfwidth1,'linear');
    
    x11=interp1(half2(1:37,1),x(1,1:37),halfwidth2,'linear');
    x22=interp1(half2(37:end,1),x(1,37:end),halfwidth2,'linear');
    
    % PLOT DECORRELATING DISTANCE
    
%      decordistALLnorm(1,:) = (decordistALL(:,1) - min(decordistALL(:,1))) / ( max(decordistALL(:,1)) - min(decordistALL(:,1)) );
%      decordistALLnorm(2,:) = (decordistALL(:,2) - min(decordistALL(:,2))) / ( max(decordistALL(:,2)) - min(decordistALL(:,2)) );

    figure(6);
    p1=plot(linspace(-120,120,73),decordistALL(:,1),'color','k','DisplayName','Control');
    hold on
%     plot([x1,x2],[halfwidth1;halfwidth1],'color','k')
    hold on
    p2=plot(linspace(-120,120,73),decordistALL(:,2),'color','r','DisplayName','PAE');
    hold on
%     plot([x11,x22],[halfwidth2;halfwidth2],'color','r')
    xlabel('Distance (cm)')
    ylabel('Correlation')
    title('Decorrelation Distance')
%     legend([p1 p2],'Control','PAE','FontSize',12)
    legend('show','Location','best')
    
    controlHalfWidth=x2-x1;
    PAEHalfWidth=x22-x11;  
end
% STATS
% BASE STATS OFF OF THESE
% Kjelstrup... Moser 2008; Battaglia, SUTHERLAND, MCNAUGHTON 2004;
% Ravassard... MEHTA 2013; Maurer...McNaughton 2005

