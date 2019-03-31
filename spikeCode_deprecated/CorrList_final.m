batchpath = '/Users/bjclark/Desktop/Dropbox/EC_mac/control_sham/BC218 - TT 5HD 2PCorBC 1GC/__data/2010-09-29_21-37-45';

f_list = FindFiles('*Correlation.txt', 'StartingDirectory', batchpath);

corrALL = [];
for i = 1:length(f_list);
        data = importdata(f_list{i});
        corrALL = [corrALL; data];
end

% NDs_append = [];
% for j = 1:length(NDs_all)
%     NDs_append = [NDs_append; NDs_all(j).numbNDs];
% end
% 


