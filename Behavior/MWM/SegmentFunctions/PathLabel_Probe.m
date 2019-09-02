%% Compile Labeled Probe Data
for timepoint=1:3
    if timepoint==3
        load('T3labelsProbe.mat'); load('T3multiProbe.mat');  %load applicable data files
        TgM=[1:7,15]; TgF=[8:14,16]; WtM=17:22; WtF=23:28;
        SEX=[1 1 1 1 1 1 1 0 0 0 0 0 0 0 1 0 1 1 1 1 1 0 0 0 0 0 0]'; %for T1 & T2
        GENOTYPE=[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0]'; %for T1 & T2
    elseif timepoint==2
        load('T2probeLabels.mat'); load('T2multiProbe.mat');  %load applicable data files
        TgM=[1:7,15]; TgF=[8:14,16]; WtM=17:22; WtF=23:28;
        SEX=[1 1 1 1 1 1 1 0 0 0 0 0 0 0 1 0 1 1 1 1 1 1 0 0 0 0 0 0]'; %for T1 & T2
        GENOTYPE=[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0]'; %for T1 & T2
    elseif timepoint==1
        load('T1probeLabels.mat'); load('T1multiProbe.mat');  %load applicable data files
        % DEFINE GROUPS
        TgM=[1:7,15]; TgF=[8:14,16]; WtM=17:22; WtF=23:28;
        SEX=[1 1 1 1 1 1 1 0 0 0 0 0 0 0 1 0 1 1 1 1 1 1 0 0 0 0 0 0]'; %for T1 & T2
        GENOTYPE=[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0]'; %for T1 & T2
    end
    
    
    % GET SEGMENT ID
    segID=[]; trialID=[]; Genotype=[]; Sex=[];
    for ii=1:length(GENOTYPE)
        Genotype=[Genotype;repmat(GENOTYPE (ii),sum(label(:,2)==ii),1)];
        Sex=[Sex; repmat(SEX(ii),sum(label(:,2)==ii),1)];
    end
    
    %COMPILE INDEX WITH LABELS
    fullLabel=[Genotype Sex label(6:end,2:3)]; %concatenate group index with data
    fullLabel=[zeros(5,4); fullLabel];%add padded zeros back to data
    
    %COMPILE MULTI STRATEGIES
    name=fieldnames(multi);
    for k=1:length(fieldnames(multi)) %cycle through fieldnames
        tempStrat=multi.(name{k});    % pull strategy data
        pullIndex=str2double(regexp(name{k,:},'(?<=i\s*)\d*', 'match'));  %pull index of strategy
        tempNewRow=fullLabel(pullIndex,:);
        for j=1:size(tempStrat,2)
            newRow(j,:)=[tempNewRow(:,1:end-1) tempStrat(:,j)];
        end
        fullLabel=[fullLabel; newRow];
    end
    
    fullLabel=sortrows(fullLabel,3);
    fullLabel(fullLabel(:,4)==12,:)=[];
    
    % COMPILE FREQ BY TRIAL AND SUBJECT FOR CONVOLUTION ANALYSIS COMPUTATIONS
    for l=1:length(unique(fullLabel(:,3)))
        temp=fullLabel(fullLabel(:,3)==l,4);
        for ii=1:11
            por(ii,1)=sum(temp(:,1)==ii)/size(temp,1);
        end
        allcollect(l,:)=por';
    end
    
    % TgM  FREQ Total
    temp=fullLabel(fullLabel(:,1)==1 & fullLabel(:,2)==1,4);
    for ii=1:11
        tgmTotal(ii,1)=sum(temp(:,1)==ii)/size(temp,1);
    end
    clear temp;
    % WtM  FREQ Total
    temp=fullLabel(fullLabel(:,1)==0 & fullLabel(:,2)==1,4);
    for ii=1:11
        wtmTotal(ii,1)=sum(temp(:,1)==ii)/size(temp,1);
    end
    clear temp;
    %  TgF FREQ Total
    temp=fullLabel(fullLabel(:,1)==1 & fullLabel(:,2)==0,4);
    for ii=1:11
        tgfTotal(ii,1)=sum(temp(:,1)==ii)/size(temp,1);
    end
    clear temp;
    %  WTF FREQ Total
    temp=fullLabel(fullLabel(:,1)==0 & fullLabel(:,2)==0,4);
    for ii=1:11
        wtfTotal(ii,1)=sum(temp(:,1)==ii)/size(temp,1);
    end
    clear temp;
    T3freqProbe=[wtmTotal tgmTotal wtfTotal tgfTotal];
    
    
    % Target-Direct
    TDFreq=sum(T3freqProbe(1:2,:),1);
    % Target-Indirect
    TIFreq=sum(T3freqProbe(3:4,:),1);
    % Spatial-Indirect
    SIFreq=sum(T3freqProbe(5:8,:),1);
    % Non-Spatial
    NSFreq=sum(T3freqProbe(9:11,:),1);
    
    CompFreq=[TDFreq; TIFreq; SIFreq; NSFreq];
    
    
    %CHI SQUARE TEST OF OVERALL STRATEGIES USED
    TgChiIdx=fullLabel(fullLabel(:,2)==1,:);
    WtChiIdx=fullLabel(fullLabel(:,2)==0,:);
    
    [TgChi,chi2,p,labels] = crosstab(TgChiIdx(:,4),TgChiIdx(:,2),TgChiIdx(:,1)) %Within Tg
    [WtChi,chi2,p,labels] = crosstab(WtChiIdx(:,4),WtChiIdx(:,2),WtChiIdx(:,1)) %Within WT
    [TgWtChi,chi2,p,labels] = crosstab(fullLabel(:,4),fullLabel(:,2)) %between Genotype
    [MFChi,chi2,p,labels] = crosstab(fullLabel(:,4),fullLabel(:,1)) %between Sex
    [TgWtMFChi,chi2,p,labels] = crosstab(fullLabel(:,4),fullLabel(:,2),fullLabel(:,1)) %between Genotype x Sex
    
    ChiPostHoc
    
    clear     segID trialID Genotype Sex GENOTYPE SEX
end

 