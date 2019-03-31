batchpath = '/Users/bjclark/Desktop/Dropbox/EC_mac/control_sham/BC225 - Cont 3HD 3BC 2GC/__data/2010-10-26_16-25-58';

f_list = dir([batchpath '/*_Correlation.mat']);

for i = 1:length(f_list);
        NDs_all(i) = load(f_list(i).name,'numbNDs');
end

NDs_append = [];
for j = 1:length(NDs_all)
    NDs_append = [NDs_append; NDs_all(j).numbNDs];
end



