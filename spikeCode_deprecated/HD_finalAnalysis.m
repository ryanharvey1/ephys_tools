batchpath = '/Users/bjclark/Desktop/Dropbox/EC_analysis/__RawData_&_Summary/__lesion_results_PaS_sm';

f_list = dir([batchpath '/*_HD_Properties.mat']);

for i = 1:length(f_list);
        r(i) = load(f_list(i).name,'mean_vector_length');
        FS_corr(i) = load(f_list(i).name,'FirstSec_Corr');
        Q_corr(i) = load(f_list(i).name,'FourQuart_MeanCorr');
        dIC(i) = load(f_list(i).name,'Direct_infoContent');
        pFR(i) = load(f_list(i).name,'peak_Firing_Rate');
        pd(i) = load(f_list(i).name,'preferred_Direction');
        or(i) = load(f_list(i).name,'overall_rate');
        halfW(i) = load(f_list(i).name,'Directional_Range_HalfWidth');
%         halfW_bins(i) = load(f_list(i).name,'Directional_Range_HalfWidth_bins');
end

r_Append = [];
for j = 1:length(r)
    r_Append = [r_Append; r(j).mean_vector_length];
end

FS_Append = [];
for k = 1:length(FS_corr)
    FS_Append = [FS_Append; FS_corr(k).FirstSec_Corr];
end

Q_Append = [];
for l = 1:length(Q_corr)
    Q_Append = [Q_Append; Q_corr(l).FourQuart_MeanCorr];
end

dIC_Append = [];
for m = 1:length(dIC)
    dIC_Append = [dIC_Append; dIC(m).Direct_infoContent];
end

pFR_Append = [];
for n = 1:length(pFR)
    pFR_Append = [pFR_Append; pFR(n).peak_Firing_Rate];
end

pd_Append = [];
for q = 1:length(pd)
    pd_Append = [pd_Append; pd(q).preferred_Direction];
end

or_Append = [];
for s = 1:length(or)
    or_Append = [or_Append; or(s).overall_rate];
end

halfW_Append = [];
for t = 1:length(halfW)
    halfW_Append = [halfW_Append; halfW(t).Directional_Range_HalfWidth];
end
% 
% halfW_bins_Append = [];
% for u = 1:length(halfW_bins)
%     halfW_bins_Append = [halfW_bins_Append; halfW_bins(u).Directional_Range_HalfWidth_bins];
% end

data_Appendi = [or_Append, r_Append, FS_Append, Q_Append, dIC_Append, pFR_Append, pd_Append halfW_Append];
dataFilt_or = find(data_Appendi(:,1) <= 10 & data_Appendi(:,2) >= 0.2829);
dataFR = find(data_Appendi(:,1) <= 10.0);
data_Append = data_Appendi(dataFR,:);

% percHD = (numel(data_Append(:,3))/numel(dataFR(:,1)))*100;

% 
% rFilt = data_Append(:,2);
% fsFilt = data_Append(:,3);
% qFilt = data_Append(:,4);









