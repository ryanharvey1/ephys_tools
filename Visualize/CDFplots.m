function [ AllStats ] = CDFplots(Group1,Group2,GroupNames,VarNames,plots)
%ScatterBox Calculates parametric & non-parametric tests and plots data
%
%   Input:
%       Group1: Column vector or matrix of data from first group
%       Group2: Column vector or matrix of data from second group
%       GroupNames: Cell array of group names. i.e.  GroupNames={'Control','PAE'}
%       VarNames: Cell array of Variable names. i.e. VarNames={'InfoContent','Coherence','PeakRate'}
%       plot: 1=plot only sig differences; 2=plot everything; 3=plot nothing
%   Output:
%       AllStats: Cell Array of Wilcoxon rank sum test and Cohen's D results
%       Plot: Creates subplot of all significant comparisons
%             Scatter dots represent raw data, orange bar represents medium value of data
%
%   Notes:
%       -This function requires notBoxPlot.m by Rob Campbell
%       -This function needs cohend.m* by Guillaume Rousselet - University of Glasgow
%       cohend.m can be found here: https://garstats.wordpress.com *(you can replace cohend.m with any cohen's d calculation)
%       -Group1 & Group2 can have different lengths
%
% Ryan E Harvey 3/29/2017

com=which('CDFplots');
com=strsplit(com,filesep);
basedir=[com{1},filesep,'Users',filesep,com{3},filesep,'GoogleDrive',filesep,'MatlabDir'];
addpath([basedir,filesep,filesep,'MatlabStatsUofG'],...
    [basedir,filesep,'RC_notBoxPlot'],...
    [basedir,filesep,'CircStat2012a'],...
    [basedir,filesep,'BClarkToolbox',filesep,'Analysis']);

% MAKE MATRIX = LENGTH WITH NANS
if size(Group1,1)~=size(Group2,1)
    if size(Group1,1)-size(Group2,1)>size(Group2,1)-size(Group1,1)
        Group2=[Group2;nan(size(Group1,1)-size(Group2,1),size(Group1,2))];
    else
        Group1=[Group1;nan(size(Group2,1)-size(Group1,1),size(Group2,2))];
    end
end

% REMOVE UNDERSCORES FROM VAR NAMES
VarNames=regexprep(VarNames,'_','','emptymatch');
if ischar(VarNames)
    VarNames={VarNames};
end

% Stats Decision Tree * WORKING *
for vars=1:size(VarNames,2)
    % Remove NaNs
    con=Group1(:,vars); 
    con(isnan(con))=[]; 
    con(isnan(con) | isinf(con))=[];
    PAE=Group2(:,vars); 
    PAE(isnan(PAE) | isinf(PAE))=[];
    % Kolmogorov-Smirnov test to test if the sample was derived from a standard normal distribution
    Kolmogorov=sum([kstest(zscore(con)),kstest(zscore(PAE))]);

    if Kolmogorov==2
        % Wilcoxon rank sum test
        [p,hypothesis1(vars),statz] = ranksum(Group1(:,vars),Group2(:,vars));
        % median
        median1=nanmedian(Group1(:,vars)); median2=nanmedian(Group2(:,vars));
        % POWER **** (Power for a non parametric test is more tricky than this--- Need to re-work)
        %         statsCON = bootstrp(length(con),@(x)[median(x),std(x)],con);
        %         statsPAE = bootstrp(length(PAE),@(x)median(x),PAE);
        %         try
        %             Power=sampsizepwr('z',[mean(statsCON(:,1)),std(statsCON(:,2))],mean(statsPAE),[],max([length(con),length(PAE)]),'Ratio',min([length(con),length(PAE)]));
        %         catch
        %             Power=NaN;
        %         end
    elseif Kolmogorov<2
        % Check for equal variance
        variances=vartest2(con,PAE);
        if variances==1
            % Two-Sample t test for unequal variance
            [hypothesis1(vars),p,ci,statz]=ttest2(Group1(:,vars),Group2(:,vars),'Vartype','unequal');
            mean1=nanmean(Group1(:,vars)); mean2=nanmean(Group2(:,vars));
        elseif variances==0
            % Two-Sample t test for equal variance
            [hypothesis1(vars),p,ci,statz]=ttest2(Group1(:,vars),Group2(:,vars));
            mean1=nanmean(Group1(:,vars)); mean2=nanmean(Group2(:,vars));
        end
        % POWER
        % Bootstrap two sample distributions with equal means and std for first group and equal means for second group
        %         statsCON = bootstrp(length(con),@(x)[mean(x) std(x)],con);
        %         Power=sampsizepwr('t2',[mean(statsCON(:,1)),mean(statsCON(:,2))],mean(bootstrp(length(PAE),@(x)mean(x),PAE)),[],max([length(con),length(PAE)]),'Ratio',min([length(con),length(PAE)]));
    end
    % Cohen's D
    cod = cohend(Group1(:,vars),Group2(:,vars));
    % Concat Stats
    if Kolmogorov==2
        AllStats(vars,:)=strcat(cellstr(VarNames(vars)),': Z=',num2str(round(statz.zval,2)),', p=',num2str(round(p,3)),', d=',num2str(round(cod,2)),...
            ',','',GroupNames(1),' Median=',num2str(round(median1,2)),', ',GroupNames(2),' Median=',num2str(round(median2,2)));
    else
        AllStats(vars,:)=strcat(cellstr(VarNames(vars)),': t(',num2str(round(statz.df)),')','=',num2str(round(statz.tstat,2)),', p=',num2str(round(p,3)),', CI[',[num2str(round(ci(1),2)),',',num2str(round(ci(2),2))],']',...
            ', d=',num2str(round(cod,2)),', ',GroupNames(1),' Mean=',num2str(round(mean1,2)),', ',GroupNames(2),' Mean=',num2str(round(mean2,2)));
    end
end
clear vars
if length(VarNames)>4; font=10;else; font=20;end
if plots~=3
    var4plot=1;
%     fig1=figure;fig1.Color=[1 1 1];
    for vars=1:size(VarNames,2)
        if hypothesis1(vars)==1 && plots==1 || plots==2
            if plots==1;SQR=round(sqrt(sum(hypothesis1))); elseif plots==2; SQR=round(sqrt(size(VarNames,2))); end
            if SQR*SQR~=size(VarNames,2)
                subplot(SQR,SQR+1,var4plot)
            else
                subplot(SQR,SQR,var4plot)
            end
            fig1.Color=[1 1 1];
%             h = notBoxPlot([Group1(:,vars),Group2(:,vars)],[],'jitter',0.6,'style','sdline','markMedian',true);
            
            [f1,x1] = ecdf(Group1(:,vars));
            [f2,x2] = ecdf(Group2(:,vars));
            
            p1=plot(x1,f1);
            set(p1,'LineWidth',4,'Color','k')
            hold on
            p2=plot(x2,f2);
            set(p2,'LineWidth',4,'Color','r')
            ylabel('Cumulative Frequency')
            xlabel(VarNames(vars))
            ax=gca;
            set(ax,'FontSize',font,'FontWeight','bold','LineWidth',2,'box','off')
            var4plot=var4plot+1;
        end
    end
end
end

