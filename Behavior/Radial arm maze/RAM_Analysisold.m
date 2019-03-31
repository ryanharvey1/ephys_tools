% RAM_Analysis

% This scrip imports data and creates stats for the 8-arm radial maze
% Currently supports 1 or 2 groups / 4 trials per session

% Input: Excel sheet with the following format

% GROUP 1                                       GROUP 2
                        % WMC WMI RMC RMI Total         % WMC WMI RMC RMI Total 
%Animal 1       % Trial 1                       % Trial 1
                % Trial 2                       % Trial 2
                % Trial 3                       % Trial 3
                % Trial 4                       % Trial 4
%Animal 2       % Trial 1                       % Trial 1
                % Trial 2                       % Trial 2
                % Trial 3                       % Trial 3
                % Trial 4                       % Trial 4
                % etc...
                
% Output: A matrix in the following format
        % Animal 1 PercentCorrect WMC WMI RMC RMI    
        % Animal 2 PercentCorrect WMC WMI RMC RMI 
        % etc...
        % &
        % ttest results for all measures
% Ryan Harvey 1/12/17        
        
clc,clear
%###########################################
% INPUT # OF ANIMALS AND RAW DATA
NumofAnimals=10;
NumofGroups=2;
data=xlsread('/Users/RyanHarvey/Dropbox/school work/UNM/Lab/Projects/Anterior Thalamic inactivation/Data/test 4:5/Test 4&5 v2',15); 
%###########################################

% CALCULATE STATS FOR GROUP #1
ianimals=1;
ianimals2=1;
session=4:4:NumofAnimals*4;
ii=1;
for i=1:length(data)
    PercentCorrect(ii,:)=(1-(sum(data(i,1:4))/data(i,5)))*100;
    WMC(ii,:)=data(i,1);
    WMI(ii,:)=data(i,2);
    RMC(ii,:)=data(i,3);
    RMI(ii,:)=data(i,4);
    ii=ii+1;
    if i==session(ianimals2)
        PercentCorrect(isnan(PercentCorrect))=0; % REMOVE NAN
        PercentCorrectKeep(ianimals2,:)=nanmean(PercentCorrect); % MEAN PERCENT
        WMCkeep(ianimals2,:)=nanmean(WMC);
        WMIkeep(ianimals2,:)=nanmean(WMI);
        RMCkeep(ianimals2,:)=nanmean(RMC);
        RMIkeep(ianimals2,:)=nanmean(RMI);
        ianimals2=ianimals2+1;
        ii=1;
        clear PercentCorrect WMC WMI RMC RMI
    end
end
Group1Stats=[PercentCorrectKeep WMCkeep WMIkeep RMCkeep RMIkeep];

if NumofGroups>1
    % CALCULATE STATS FOR GROUP #2
    ianimals=1;
    ianimals2=1;
    session=4:4:NumofAnimals*4;
    ii=1;
    for i=1:length(data)
        PercentCorrect(ii,:)=(1-(sum(data(i,6:10))/data(i,10)))*100;
        WMC(ii,:)=data(i,6);
        WMI(ii,:)=data(i,7);
        RMC(ii,:)=data(i,8);
        RMI(ii,:)=data(i,9);
        ii=ii+1;
        if i==session(ianimals2)
            PercentCorrect(isnan(PercentCorrect))=0; % REMOVE NAN
            PercentCorrectKeep(ianimals2,:)=nanmean(PercentCorrect); % MEAN PERCENT
            WMCkeep(ianimals2,:)=nanmean(WMC);
            WMIkeep(ianimals2,:)=nanmean(WMI);
            RMCkeep(ianimals2,:)=nanmean(RMC);
            RMIkeep(ianimals2,:)=nanmean(RMI);
            ianimals2=ianimals2+1;
            ii=1;
            clear PercentCorrect WMC WMI RMC RMI
        end
    end
    Group2Stats=[PercentCorrectKeep WMCkeep WMIkeep RMCkeep RMIkeep];
    
    variables={'PercentCorrect' 'WMC' 'WMI' 'RMC' 'RMI'};
    for i=1:length(variables)
        [h,p,ci,stats] = ttest(Group1Stats(:,i),Group2Stats(:,i));   
        
        d=mean((Group2Stats(:,i))-mean(Group1Stats(:,i)))/sqrt((std(Group1Stats(:,i))^2+std(Group2Stats(:,i))^2)/2); 
        
        staaaats(i,:)=strcat(cellstr(variables(i)),': t(',num2str(stats.df),...
            ')=',num2str(stats.tstat),', p=',(num2str(p)),'  Confidence: [',...
            num2str(ci(1)),',',num2str(ci(2)),']','  sd: ', num2str(stats.sd),'  Cohen''s d: ',num2str(d));
    end
end
staaaats



