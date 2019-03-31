%%% code for second round of revisions for paper

batchpath = '/Users/bjclark/Desktop/Dropbox/EC_analysis/__RawData_&_Summary/__RateMapData/__lesion_data_26x26/sm_mec';
% batchpath2 = '/Users/bjclark/Desktop/Dropbox/EC_analysis/__RawData_&_Summary/__control_results_files_MEC';

f_list = dir([batchpath '/*_RateMapStats.mat']);
h_list = dir([batchpath '/*_HD_Properties.mat']);

for i = 1:length(f_list);
        GS1(i) = load(f_list(i).name,'maxGS1');
        SinuGS(i) = load(f_list(i).name,'maxSinuGrid');
        IC(i) = load(f_list(i).name,'InformationContent');
        Spars(i) = load(f_list(i).name,'Sparsity');
        pFR(i) = load(f_list(i).name,'peakFR');
%         GS2(i) = load(f_list(i).name,'maxGS2');
%         GS3(i) = load(f_list(i).name,'maxGS3');
end

for h = 1:length(h_list);
        ovFR(h) = load(h_list(h).name,'overall_rate');
end

GS1_Append = [];
for j = 1:length(GS1)
    GS1_Append = [GS1_Append; GS1(j).maxGS1];
end

SinuGS_Append = [];
for k = 1:length(SinuGS)
    SinuGS_Append = [SinuGS_Append; SinuGS(k).maxSinuGrid];
end

IC_Append = [];
for m = 1:length(IC)
    IC_Append = [IC_Append; IC(m).InformationContent];
end

Spars_Append = [];
for n = 1:length(Spars)
    Spars_Append = [Spars_Append; Spars(n).Sparsity];
end

pFR_Append = [];
for q = 1:length(pFR)
    pFR_Append = [pFR_Append; pFR(q).peakFR];
end
% 
% mSGS_Append = [];
% for r = 1:length(mSinuGS)
%     mSGS_Append = [mSGS_Append; mSinuGS(r).meanSinuGrid];
% end

ovFR_Append = [];
for s = 1:length(ovFR)
    ovFR_Append = [ovFR_Append; ovFR(s).overall_rate];
end

% data = [ovFR_Append GS1_Append SinuGS_Append IC_Append Spars_Append pFR_Append];
% datai = find(data(:,1) <= 10.0 & data(:,3) >= 0.4393); 
% dataFR = find(data(:,1) <= 10.0);
% dataGSFILT = data(dataFR,:);
% 
% % percGS1 = (numel(dataGSFILT(:,3))/numel(dataFR(:,1)))*100;
% % 
% e = -1:0.05:1.5;
% [n,xout] = hist(dataGSFILT(:,3),e);
% freq = ((n/sum(n))*100);
% CumFreqGS1 = cumsum(hist(dataGSFILT(:,3),50));

% 
% bar(xout,freq,0.8,'k');
% % 
% plot(xout,freq,'LineWidth',1.5,'color','r');
% ylim([0 15]);
% xlim([-0.5 1.0]);
% box off
% 
% scatter(dataGSFILT(:,7),dataGSFILT(:,3),'filled','r');
% ylim([-0.5 1.5]);





