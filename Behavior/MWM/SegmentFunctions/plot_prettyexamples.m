

load('C:\Users\Ben Clark''s Lab\Google Drive\MATLAB\MWM\Stats\MovementTracking\T1FreqCompiled.mat')
keep fullLabel
load('C:\Users\Ben Clark''s Lab\Google Drive\MATLAB\MWM\Stats\PathFeatures\pathResults_092018.mat')

for i=1:unique(fullLabel(:,4))
    figure; 
    labelTemp=sortrows(fullLabel(fullLabel(:,4)==i,:),5);
    for ii=1:length(labelTemp)
        plot(pathResults.T1_t1_sbj1.segments(ii).segments(:,1),pathResults.T1_t1_sbj1.segments(ii).segments(:,2)); hold on;
    end
end