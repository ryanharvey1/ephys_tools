function [table,SumofNaN1,SumofNaN2,PercentofNaN] = NonDetects(table)
%NonDetects Summary of this function goes here
%   Detailed explanation goes here

 % IF TRACKING NON-DETECT>>>REPLACE COORDINATE WITH NAN
    table(table==0)=NaN;
    % CALC NUMBER AND PERCENT OF NON-DETECTS
    SumofNaN1=sum(isnan(table(:,2))); SumofNaN2=sum(isnan(table(:,4)));
    PercentofNaN=(max([SumofNaN1 SumofNaN2])/length(table))*100;
    disp(['Filling in ', num2str(max([SumofNaN1 SumofNaN2])), ' non-detects, which is ', num2str(PercentofNaN), ' percent of your data.'])
    
    % IF *FIRST* COORDINATE IS NAN, USE AVERAGE OF FIRST COORDINATE
    for icolumn=2:5
        if isnan(table(1,icolumn))==1
            for inan=1:length(table)
                if isnan(table(inan+1,icolumn))==0 % COUNTS # OF ROWS OF NANs
                    n=inan+1; % # OF NANS
                    break
                end
            end
            table(1,icolumn)=nanmean(table(1:n,icolumn)); % AVG OF n
        end
        % IF *LAST* COORDINATE IS NAN, USE AVERAGE OF LAST COORDINATE
        if isnan(table(end,icolumn))==1
            for inan=length(table):-1:1
                if isnan(table(inan-1,icolumn))==0 % COUNTS # OF ROWS OF NANs
                    n=inan-1; % # OF NANS
                    break
                end
            end
            table(end,icolumn)=nanmean(table(n:end,icolumn)); % AVG OF n
        end
        [table]=FillNaN2(table,icolumn);
    end
end

