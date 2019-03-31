function [PhaseStats] = PhaseLock(Group1Deg,Group2Deg,Group1RayP,Group2RayP)
%PhaseLock Calculates Phase locking differences between two groups
%   
%   Input: 
%         Group1Deg: Group 1 number of cells (x rows) BY number of frequencies (6 columns) matrix of mean degrees
%         Group2Deg: Group 2 number of cells (x rows) BY number of frequencies (6 columns) matrix of mean degrees
%         Group1RayP: Group 1 number of cells (x rows) BY number of frequencies (6 columns) matrix of Rayleigh PValues
%         Group2RayP: Group 2 number of cells (x rows) BY number of frequencies (6 columns) matrix of Rayleigh PValues
%
%   Output:
%         PhaseStats: Non parametric multi-sample test for equal medians. Similar to a
%                     Kruskal-Wallis test for linear data.
%                     Watson's statistic for the Goodness-of-fit tests on a circle  
%         Plots: Signigicanly Phase locked mean degrees on Phase
%
%   Example:
%         [PhaseStats] = PhaseLock(control(:,(41:46)),PAE(:,(41:46)),control(:,(29:34)),PAE(:,(29:34)));
%
% Ryan E Harvey 3/29/2017


variables={'DELTA','THETA','ALPHA','BETA','GAMMA','HI-GAMMA'};
% ANOVAStats=[];
for freq=1:6
    % Filter by sig phase locking
%     Group1DegWork=Group1Deg(Group1RayP(:,freq)<=0.05,freq); Group2DegWork=Group2Deg(Group2RayP(:,freq)<=0.05,freq);
    Group1DegWork=Group1Deg(:,freq); Group2DegWork=Group2Deg(:,freq);

    % Calculate Stats
%     % Watson-Williams
%     [~,ANOVATab] = circ_wwtest(deg2rad(Group1DegWork(:,1)),deg2rad(Group2DegWork(:,1)));
%     ANOVAStats=[ANOVAStats;[variables(freq);'-';'-';'-'],ANOVATab];
%     
    % Watson's Goodness-of-fit tests 
    [U2,p] = watson1962(deg2rad(Group1DegWork(:,1)),deg2rad(Group2DegWork(:,1)));
    
    % Kruskal-Wallis
    [pval,med,P] = circ_cmtest(deg2rad(Group1DegWork(:,1)),deg2rad(Group2DegWork(:,1)));

    PhaseStats(freq,:)=[strcat(cellstr(variables(freq)),': ','Kruskal-Wallis Statistic:',num2str(P),'  Shared Medium:',num2str(rad2deg(med)),'  PVal:',num2str(pval)),...
        strcat('Watson''s Goodness of fit: ','  U2:',num2str(U2),'  Pvalue:',num2str(p))];
    clear pval med P p U2
end

subNum=1;
for iii=1:2
    if iii==1;group=Group1Deg; elseif iii==2;group=Group2Deg; end;
    ii=1;
    for i=1:6
%         if iii==1; GroupWorking=group(Group1RayP(:,i)<=0.05,i); end
%         if iii==2; GroupWorking=group(Group2RayP(:,i)<=0.05,i); end
        if iii==1; GroupWorking=group(:,i); end
        if iii==2; GroupWorking=group(:,i); end
        figure(2); subplot(2,6,subNum)
        % plot wave
        steps=(2*pi)/(length(GroupWorking)-1);
        t = 0:steps:(2*pi);
        y = sin(t);
        hh=plot(t,y,'k'); hold on
        hh.LineWidth=5;
        set(gca, 'XTick', []);
        set(gca, 'YTick', []);
        box off
        title(variables(ii),'FontSize',12); ii=ii+1;
        % scatter spikes
        spikes=sin(deg2rad(GroupWorking(:,1))');
        h=scatter(deg2rad(GroupWorking(:,1))',spikes,'Marker','x');
        h.MarkerEdgeColor=[1 0 0];
        h.MarkerFaceColor=[1 0 0];
        h.LineWidth = 8;
        axis([0 2*pi -1.2 1.2])
        hold on
        % Plot mean Spike
        meanspike=circ_mean(deg2rad(GroupWorking));
        if meanspike<0; meanspike=deg2rad(rad2deg(meanspike)+360); end
        spikes=sin(meanspike);
        h2=scatter(meanspike,spikes,'Marker','o');
        hold on
        h2.MarkerEdgeColor='b';
        h2.MarkerFaceColor='b';
        h2.LineWidth = 10;
        subNum=subNum+1;
    end
end
end

