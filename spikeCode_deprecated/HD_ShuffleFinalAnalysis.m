batchpath = '/Users/bjclark/Desktop/Dropbox/EC_analysis/__RawData_&_Summary/__control_results_files_only';

f_list = dir([batchpath '/*_HD_shuff.mat']);

for i = 1:length(f_list);
        rSHUFFLEi(i) = load(f_list(i).name,'rSHUFF');
        corrSHUFFLEi(i) = load(f_list(i).name,'corrSHUFF');
end

rSHUFFLE = [];
for j = 1:length(rSHUFFLEi)
    rSHUFFLE = [rSHUFFLE; rSHUFFLEi(j).rSHUFF];
end

corrSHUFFLE = [];
for k = 1:length(corrSHUFFLEi)
    corrSHUFFLE = [corrSHUFFLE; corrSHUFFLEi(k).corrSHUFF];
end

percNinetyFive_MVR = prctile(rSHUFFLE,95);
percNinetyNine_MVR = prctile(rSHUFFLE,99);
figure (1), hist(rSHUFFLE,100);
hold on
% line([percNinetyNine_MVR percNinetyNine_MVR], [0 18000], 'LineWidth',1, 'color', 'r')
line([percNinetyFive_MVR percNinetyFive_MVR], [0 30000], 'LineWidth',1.5, 'color', 'r')
title('Shuffled Distribution of Mean Vector Length (red = 95th)', 'fontSize',18);
box off;

% shuffled ditribution for corr
percNinetyFive_CORR = prctile(corrSHUFFLE,95);
percNinetyNine_CORR = prctile(corrSHUFFLE,99);
figure (2), hist(corrSHUFFLE,100);
hold on
% line([percNinetyNine_CORR percNinetyNine_CORR], [0 15000], 'LineWidth',1, 'color', 'r')
line([percNinetyFive_CORR percNinetyFive_CORR], [0 25000], 'LineWidth',1.5, 'color', 'r')
title('Shuffled Distribution of Mean Corr (red = 95th)', 'fontSize',18);
box off;







