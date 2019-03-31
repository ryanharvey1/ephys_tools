% numberofspikes
% Finally, different smoothing parameters were used for cells with different 
% numbers of spikes (smoothing filter passed 5x in one case and 10x in another). 
% Again, where there differences in the numbers of cells that fell into the 
% low vs higher smoothing categories across the two groups: tilted vs control? 
% If so, some evidence would need to be provided that such differences didn't 
% actually cause the reported effects.

% - CALCULATE NUMBER OF SPIKES IN EACH GROUP AND SHOW THAT MOST ARE ABOVE 100 SPIKES SO DIFFERENT SMOOTHING PARAMETERS WOULD HAVE NO EFFECT


load('/Users/ryanharvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/NewOCC_Map_workspace2.mat')
    
DispVarNames={'PeakRate','nSpikes','OverallFR','NumbActiveBins','sparsity','InformationContent','Coherence','Field2Wall','borderScore','FieldWidth','infieldFR','outfieldFR','E','c','p'};

%  DIFFERENCES IN SPIKES BETWEEN GROUPS AND MIN NUMBERS
for i=1:5
    ScatterBox(resultsC(:,2,i),resultsT(:,2,i),{'Con','Con'},DispVarNames(2),2)
    disp(['Session ',num2str(i),'Control Min: ',num2str(min(resultsC(:,2,i))),' Tilted Min: ',num2str(min(resultsT(:,2,i)))])
end

% PERCENT BELOW 100
for i=1:5
    [h,p, chi2stat,df] = prop_test([sum(resultsC(:,2,i)<=100), sum(resultsT(:,2,i)<=100)],[length(resultsC(:,2,i)), length(resultsT(:,2,i))], 0);

    disp(['(Session ',num2str(i),')',' Control below 100: ',num2str(sum(resultsC(:,2,i)<=100)/length(resultsC(:,2,i))),' Tilted below 100: ',num2str(sum(resultsT(:,2,i)<=100)/length(resultsT(:,2,i))),' p= ',num2str(p)])
end

