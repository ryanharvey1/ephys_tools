batchpath = '/Users/bjclark/Desktop/Dropbox/EC_analysis/__RawData_&_Summary/__lesion_results_files_PaS';

f_list = dir([batchpath '/*_NDs.mat']);

for i = 1:length(f_list);
        NDs_all(i) = load(f_list(i).name,'numbNDs');
end

NDs_append = [];
for j = 1:length(NDs_all)
    NDs_append = [NDs_append; NDs_all(j).numbNDs];
end

