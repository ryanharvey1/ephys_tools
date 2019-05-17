function AllStats  = stat_plot(Group1,Group2,GroupNames,VarNames,varargin)
%stat_plot Calculates parametric & non-parametric tests and plots data
%
%   Input:
%       Group1: Column vector or matrix of data from first group
%       Group2: Column vector or matrix of data from second group
%       GroupNames: Cell array of group names. i.e.  GroupNames={'Control','PAE'}
%       VarNames: Cell array of Variable names. i.e. VarNames={'InfoContent','Coherence','PeakRate'}
%   Options:
%       plots: 1=plot only sig differences (default); 2=plot everything; 3=plot nothing
%       plottype: 'cdf'(default) or 'beeswarm'
%   Output:
%       AllStats: Cell Array of Wilcoxon rank sum test and Cohen's D results
%       Plot: Creates subplot of all significant comparisons
%
%
%   Notes:
%       -This function requires notBoxPlot.m by Rob Campbell
%       -This function needs cohend.m* by Guillaume Rousselet - University of Glasgow
%       cohend.m can be found here: https://garstats.wordpress.com *(you can replace cohend.m with any cohen's d calculation)
%       -Group1 & Group2 can have different lengths
%
% Ryan E Harvey 2019
%
p = inputParser;
p.addParameter('plots',1);
p.addParameter('plottype','cdf');
p.parse(varargin{:});

plots = p.Results.plots;
plottype = p.Results.plottype;


Group1(isinf(Group1))=NaN;
Group2(isinf(Group2))=NaN;

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
    % Kolmogorov-Smirnov test to test if the sample was derived from a standard normal distribution
    Kolmogorov=sum([kstest(Group1(:,vars)),kstest(Group2(:,vars))]);
    
    if Kolmogorov>=1
        % Wilcoxon rank sum test
        [p,hypothesis1(vars),statz] = ranksum(Group1(:,vars),Group2(:,vars));
        % median
        median1=nanmedian(Group1(:,vars)); median2=nanmedian(Group2(:,vars));
        
    elseif Kolmogorov==0
        % Check for equal variance
        variances=vartest2(Group1(:,vars),Group2(:,vars));
        if variances==1
            % Two-Sample t test for unequal variance
            [hypothesis1(vars),p,ci,statz]=ttest2(Group1(:,vars),Group2(:,vars),'Vartype','unequal');
            mean1=nanmean(Group1(:,vars)); mean2=nanmean(Group2(:,vars));
        elseif variances==0
            % Two-Sample t test for equal variance
            [hypothesis1(vars),p,ci,statz]=ttest2(Group1(:,vars),Group2(:,vars));
            mean1=nanmean(Group1(:,vars)); mean2=nanmean(Group2(:,vars));
        end
    end
    % Cohen's D
    cod = cohend(Group1(:,vars),Group2(:,vars));
    % Concat Stats
    if Kolmogorov>=1
        AllStats(vars,:)=strcat(cellstr(VarNames(vars)),': Z=',num2str(round(statz.zval,2)),', p=',num2str(round(p,3)),', d=',num2str(round(cod,2)),...
            ',','',GroupNames(1),' Median=',num2str(round(median1,2)),', ',GroupNames(2),' Median=',num2str(round(median2,2)));
    else
        AllStats(vars,:)=strcat(cellstr(VarNames(vars)),': t(',num2str(round(statz.df)),')','=',num2str(round(statz.tstat,2)),', p=',num2str(round(p,3)),', CI[',[num2str(round(ci(1),2)),',',num2str(round(ci(2),2))],']',...
            ', d=',num2str(round(cod,2)),', ',GroupNames(1),' Mean=',num2str(round(mean1,2)),', ',GroupNames(2),' Mean=',num2str(round(mean2,2)));
    end
end

% Plotting
clear vars
if length(VarNames)>4; font=10;else; font=20;end
if plots~=3
    var4plot=1;
    for vars=1:size(VarNames,2)
        if hypothesis1(vars)==1 && plots==1 || plots==2
            if plots==1;SQR=round(sqrt(sum(hypothesis1))); elseif plots==2; SQR=round(sqrt(size(VarNames,2))); end
            if SQR*SQR~=size(VarNames,2)
                subplot(SQR,SQR+1,var4plot)
            else
                subplot(SQR,SQR,var4plot)
            end
            if strcmp(plottype,'cdf')
                [f1,x1] = ecdf(Group1(:,vars));
                [f2,x2] = ecdf(Group2(:,vars));
                
                p1=plot(x1,f1);
                set(p1,'LineWidth',4,'Color','k')
                hold on
                p2=plot(x2,f2);
                set(p2,'LineWidth',4,'Color','r')
                ylabel('Cumulative Frequency')
                xlabel(VarNames(vars))
            elseif strcmp(plottype,'beeswarm')
                plotspread_wrapper(Group1(:,vars),Group2(:,vars),GroupNames)
                title(VarNames(vars))
            end
            ax=gca;
            set(ax,'FontSize',font,'FontWeight','bold','LineWidth',2,'box','off')
            var4plot=var4plot+1;
        end
    end
end
end

