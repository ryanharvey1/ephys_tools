clear

path='D:\Projects\PAE_PlaceCell\ProcessedData\';

sessions=dir([path,'*.mat']);
sessions={sessions.name}';


path='D:\Projects\PAE_PlaceCell\AnimalMetadata\';
rats=dir([path,'*.mat']);
rats={rats.name}';
sessions_w_tet=[];
for i=1:length(rats)
    load(rats{i})
    disp(rats{i})
    
    current_sess=sessions(contains(sessions,extractBefore(rats{i},'_')));
    nch=AnimalMetadata.ExtracellEphys.Probes.nchannels;

    for s=1:length(current_sess)
        session=extractBetween(current_sess,'_','.mat');
        area=AnimalMetadata.RecordingLogs.(session{s}).RecordingArea;
        if ischar(area) || isempty(area)
            sessions_w_tet=[sessions_w_tet;[repmat(current_sess(s),nch,1),num2cell(1:nch)',repmat({area},nch,1)]];
        else
            sessions_w_tet=[sessions_w_tet;[repmat(current_sess(s),nch,1),num2cell(1:nch)',area']];
        end
    end
    

end
combined_ID=sessions_w_tet;

% data=compileResults('D:\Projects\PAE_PlaceCell\ProcessedData');
% 
control={'RH13','RH14','LS21','LS23','LE2821','LE2823','LEM3116','LEM3120'};
pae={'RH11','RH16','LS17','LS19','LE2813','LEM3124'};
% 
% %% COMPILE GROUPS
% data.control.id=[];
% for i=1:length(control)
%     data.control.id=cat(1,data.control.id,data.(control{i}).id);
% end
% 
% data.pae.id=[];
% for i=1:length(pae)
%     data.pae.id=cat(1,data.pae.id,data.(pae{i}).id);
% end
% 
% %% COMPILE GROUP IDS
% group1id=data.control.id;
% group2id=data.pae.id;
% 
% %GET RECORDING AREA FOR ALL SESSIONS (KUBIE) and/or TETRODES (HALO)
% group1id=get_region_id(group1id,'D:\Projects\PAE_PlaceCell\AnimalMetadata');
% group2id=get_region_id(group2id,'D:\Projects\PAE_PlaceCell\AnimalMetadata');

% combined_ID=[group1id; group2id];

%%
cd D:\Projects\PAE_PlaceCell\CFCResults
path='D:\Projects\PAE_PlaceCell\CFCResults\';
% rats=dir('*.mat');
% rats={rats.name};

%INITIALIZE DATA STRUCTURES

% tempResults.controls.sess{1}=[];
% tempResults.controls.sess{2}=[];
% tempResults.controls.sess{3}=[];
% tempResults.controls.sess{4}=[];
% 
% tempResults.pae.sess{1}=[];
% tempResults.pae.sess{2}=[];
% tempResults.pae.sess{3}=[];
% tempResults.pae.sess{4}=[];
DataMat_Ses1=table;
DataMat_Ses2=table;
DataMat_Ses3=table;

rats=dir([path,'*.mat']);
rats={rats.name}';

for i=1:length(rats)
    
    cfc=load(rats{i});

    for sess=1:size(cfc.cfc,2)
        
        peakTheta=cfc.cfc{1, sess}.theta_ch.f_peak;
        
        thetaChannel=cfc.cfc{1, sess}.theta_ch.active_lfp_ind;
        
       
        area=combined_ID{ismember(combined_ID(:,1),rats{i,1}) &...
            thetaChannel==[combined_ID{:,2}]',3};
        
        areaIdx={area};
        
        ratID=rats(i);
       
        LG_MI=cfc.cfc{1, sess}.LGtortCFC.modindex;
        LG_MI_thres=cfc.cfc{1, sess}.LGtortCFC.MIthres;
        LG_MI_sig=LG_MI>LG_MI_thres;
        
        
        HG_MI=cfc.cfc{1, sess}.HGtortCFC.modindex;
        HG_MI_thres=cfc.cfc{1, sess}.HGtortCFC.MIthres;
        HG_MI_sig=HG_MI>HG_MI_thres;
        
        %Build Data Structure for whole session MI
        
        if sess==1
            DataMat_Ses1.ratID{i}= ratID;
            DataMat_Ses1.peakTheta{i}=  peakTheta;
            DataMat_Ses1.thetaChannel{i}=  thetaChannel;
            DataMat_Ses1.areaIdx{i}= areaIdx ;
            DataMat_Ses1.LG_MI{i}=  LG_MI;
            DataMat_Ses1.LG_MI_thres{i}=  LG_MI_thres;
            DataMat_Ses1.LG_MI_sig{i}=  LG_MI_sig;
            DataMat_Ses1.HG_MI{i} = HG_MI;
            DataMat_Ses1.HG_MI_thres{i}=  HG_MI_thres;
            DataMat_Ses1.HG_MI_sig{i} =HG_MI_sig;
     
        elseif sess==2
            DataMat_Ses2.ratID{i}= ratID;
            DataMat_Ses2.peakTheta{i}=  peakTheta;
            DataMat_Ses2.thetaChannel{i}=  thetaChannel;
            DataMat_Ses2.areaIdx{i}= areaIdx ;
            DataMat_Ses2.LG_MI{i}=  LG_MI;
            DataMat_Ses2.LG_MI_thres{i}=  LG_MI_thres;
            DataMat_Ses2.LG_MI_sig{i}=  LG_MI_sig;
            DataMat_Ses2.HG_MI{i} = HG_MI;
            DataMat_Ses2.HG_MI_thres{i}=  HG_MI_thres;
            DataMat_Ses2.HG_MI_sig{i} =HG_MI_sig;
        elseif sess==3
            DataMat_Ses3.ratID{i}= ratID;
            DataMat_Ses3.peakTheta{i}=  peakTheta;
            DataMat_Ses3.thetaChannel{i}=  thetaChannel;
            DataMat_Ses3.areaIdx{i}= areaIdx ;
            DataMat_Ses3.LG_MI{i}=  LG_MI;
            DataMat_Ses3.LG_MI_thres{i}=  LG_MI_thres;
            DataMat_Ses3.LG_MI_sig{i}=  LG_MI_sig;
            DataMat_Ses3.HG_MI{i} = HG_MI;
            DataMat_Ses3.HG_MI_thres{i}=  HG_MI_thres;
            DataMat_Ses3.HG_MI_sig{i} =HG_MI_sig;
        end
        
    end
    clear cfc
end

%% Delete empties from datamat
DataMat_Ses1(all(cellfun(@isempty,DataMat_Ses1{:,:}),2),:)=[];
DataMat_Ses2(all(cellfun(@isempty,DataMat_Ses2{:,:}),2),:)=[];
DataMat_Ses3(all(cellfun(@isempty,DataMat_Ses3{:,:}),2),:)=[];



clear sess i

%% WHOLE SESSION MI ANALYSIS


for sess=1:3
    disp(['Session ', num2str(sess), ' Analysis'])
    VarNames={'Peak Theta Frequency','MI- Theta/Low Gamma','MI- Theta/Low Gamma'};
    GroupNames={'Sacc','PAE'};
    
    
    if sess==1
        Group1_CA1=DataMat_Ses1(contains([DataMat_Ses1.ratID{:}],control) & contains([DataMat_Ses1.areaIdx{:}],'ca1'),:);
        Group2_CA1=DataMat_Ses1(contains([DataMat_Ses1.ratID{:}],pae) & contains([DataMat_Ses1.areaIdx{:}],'ca1'),:);
        
        Group1_CA3=DataMat_Ses1(contains([DataMat_Ses1.ratID{:}],control) & contains([DataMat_Ses1.areaIdx{:}],'ca3') ,:);
        Group2_CA3=DataMat_Ses1(contains([DataMat_Ses1.ratID{:}],pae) & contains([DataMat_Ses1.areaIdx{:}],'ca3'),:);
    elseif sess==2
        Group1_CA1=DataMat_Ses2(contains([DataMat_Ses2.ratID{:}],control) & contains([DataMat_Ses2.areaIdx{:}],'ca1'),:);
        Group2_CA1=DataMat_Ses2(contains([DataMat_Ses2.ratID{:}],pae) & contains([DataMat_Ses2.areaIdx{:}],'ca1'),:);
        
        Group1_CA3=DataMat_Ses2(contains([DataMat_Ses2.ratID{:}],control) & contains([DataMat_Ses2.areaIdx{:}],'ca3') ,:);
        Group2_CA3=DataMat_Ses2(contains([DataMat_Ses2.ratID{:}],pae) & contains([DataMat_Ses2.areaIdx{:}],'ca3'),:);
    elseif sess==3
        Group1_CA1=DataMat_Ses3(contains([DataMat_Ses3.ratID{:}],control) & contains([DataMat_Ses3.areaIdx{:}],'ca1'),:);
        Group2_CA1=DataMat_Ses3(contains([DataMat_Ses3.ratID{:}],pae) & contains([DataMat_Ses3.areaIdx{:}],'ca1'),:);
        
        Group1_CA3=DataMat_Ses3(contains([DataMat_Ses3.ratID{:}],control) & contains([DataMat_Ses3.areaIdx{:}],'ca3') ,:);
        Group2_CA3=DataMat_Ses3(contains([DataMat_Ses3.ratID{:}],pae) & contains([DataMat_Ses3.areaIdx{:}],'ca3'),:);
    end
    
    [group1ca1,~]=extract_double_from_table(Group1_CA1);
    [group2ca1,~]=extract_double_from_table(Group2_CA1);
    [group1ca3,~]=extract_double_from_table(Group1_CA3);
    [group2ca3,varnames]=extract_double_from_table(Group2_CA3);

%     %low gamma
%     LG_propSig_group1=sum(group1ca1(:,5))/size(group1ca1(:,5),1);
%     LG_propSig_group2=sum(group2ca1(:,5))/size(group2ca1(:,5),1);
    
    %test proportion of sessions that exhibit phase-amplitude coupling
    disp('CA1 proportion of low gamma phase-amp coupling')
    [~,p, chi2stat,df] = prop_test([sum(group1ca1(:,5))/size(group1ca1(:,5),1) sum(group2ca1(:,5))/size(group2ca1(:,5),1)] , [1 1], 0);
    disp(['X(',num2str(df),')=',num2str(chi2stat),', p=',num2str(p)])
    
    disp('CA1 proportion of high gamma phase-amp coupling')
    [~,p, chi2stat,df] = prop_test([sum(group1ca1(:,8))/size(group1ca1(:,8),1) sum(group2ca1(:,8))/size(group2ca1(:,8),1)] , [1 1], 0);
    disp(['X(',num2str(df),')=',num2str(chi2stat),', p=',num2str(p)])
    
    disp('CA3 proportion of low gamma phase-amp coupling')
    [~,p, chi2stat,df] = prop_test([sum(group1ca3(:,5))/size(group1ca3(:,5),1) sum(group2ca3(:,5))/size(group2ca3(:,5),1)] , [1 1], 0);
    disp(['X(',num2str(df),')=',num2str(chi2stat),', p=',num2str(p)])
    
    disp('CA3 proportion of high gamma phase-amp coupling')
    [~,p, chi2stat,df] = prop_test([sum(group1ca3(:,8))/size(group1ca3(:,8),1) sum(group2ca3(:,8))/size(group2ca3(:,8),1)] , [1 1], 0);
    disp(['X(',num2str(df),')=',num2str(chi2stat),', p=',num2str(p)])
    
    
%     disp(['Results for proportion of sessions that contain significant PAC',area,' session ',num2str(sess)])
%     if p<.05
%         disp(['The total number of sessions that exhibit significant'...
%             'phase-amplitude coupling is different between groups, X(',num2str(df),')=',num2str(chi2stat),', p=',num2str(p)])
%     else
%         disp(['Both groups exhibit a similar proportion of sessions with signficant '...
%             'phase-amplitude coupling, X(',num2str(df),')=',num2str(chi2stat),', p=',num2str(p)])
%     end
    
    %test magnitude difference of MI between groups
%     LG_group1=group1ca1(group1ca1(:,5)==1,3);
%     LG_group2=group2ca1(group2ca1(:,5)==1,3);
    
    %test variance
%     if ~isempty(LG_group2) && ~isempty(LG_group1)
        
%         h = vartest2(LG_group1,LG_group2);
        figure;
        LGModulationIndex_ca1=CDFplots(group1ca1(group1ca1(:,5)==1,3),...
            group2ca1(group2ca1(:,5)==1,3),{'Sacc','PAE'},{'LG Modulation Index CA1'},2);
        title(extractBetween(LGModulationIndex_ca1,',',', CI'))
                figure;

        LGModulationIndex_ca3=CDFplots(group1ca3(group1ca3(:,5)==1,3),...
            group2ca3(group2ca3(:,5)==1,3),{'Sacc','PAE'},{'LG Modulation Index CA3'},2);
        title(extractBetween(LGModulationIndex_ca3,',',', CI'))
                figure;

        HGModulationIndex_ca1=CDFplots(group1ca1(group1ca1(:,8)==1,6),...
            group2ca1(group2ca1(:,8)==1,6),{'Sacc','PAE'},{'HG Modulation Index CA1'},2);
        title(extractBetween(HGModulationIndex_ca1,',',', CI'))
                figure;

        HGModulationIndex_ca3=CDFplots(group1ca3(group1ca3(:,8)==1,6),...
            group2ca3(group2ca3(:,8)==1,6),{'Sacc','PAE'},{'HG Modulation Index CA3'},2);
        title(extractBetween(HGModulationIndex_ca3,',',', CI'))

%         if h==0 %run t-test with equal variances assumed
%             [~,p,~,LGstats] = ttest2(LG_group1,LG_group2,'Vartype','equal');
%         else %run t-test controlling for unequal variances
%             [~,p,~,LGstats] = ttest2(LG_group1,LG_group2,'Vartype','unequal');
%         end
        
%         LGmu_group1=nanmean(Group1(Group1(:,3)==1,2));
%         LGmu_group2=nanmean(Group2(Group2(:,3)==1,2));
%         
%         LGsem_group1=std(Group1(Group1(:,3)==1,2))/sqrt(size(Group1(Group1(:,3)==1,2),1));
%         LGsem_group2=std(Group2(Group2(:,3)==1,2))/sqrt(size(Group2(Group2(:,3)==1,2),1));
%         
%         disp(['Low Gamma Magnitude of MI in area ',area{region},' session ',num2str(sess)])
%         if p<.05 && LGmu_group1<LGmu_group2
%             disp(['PAE rats exhibited higher mean MI, mean=',num2str(LGmu_group2),...
%                 'sem=',num2str(LGsem_group2),' across the entire session'...
%                 'relative to Sacc controls, mean=',num2str(LGmu_group1),'sem=',num2str(LGsem_group1),'p=',num2str(p)])
%         elseif p<.05 && LGmu_group1>LGmu_group2
%             disp(['Sacc rats exhibited higher mean MI, mean=',num2str(LGmu_group1),...
%                 'sem=',num2str(LGsem_group1),' across the entire session'...
%                 'relative to PAE rats, mean=',num2str(LGmu_group2),'sem=',num2str(LGsem_group2),'p=',num2str(p)])
%         else
%             disp(['There is no significant difference between MI between groups, p=',num2str(p)])
%         end
        
%     end
    
%     
%     %high gamma
%     HG_propSig_group1=sum(Group1(:,5))/size(Group1(:,4),1);
%     HG_propSig_group2=sum(Group2(:,5))/size(Group2(:,4),1);
%     
%     %test proportion of sessions that exhibit phase-amplitude coupling
%     [~,p, chi2stat,df] = prop_test([HG_propSig_group1 HG_propSig_group2] , [1 1], 0);
%     
%     disp(['Results for proportion of sessions that contain significant PAC:',area{region},' session ',num2str(sess)])
%     if p<.05
%         disp(['The total number of sessions that exhibit significant'...
%             'phase-amplitude coupling is different between groups,X(',num2str(df),')=',num2str(chi2stat),', p=',num2str(p)])
%     else
%         disp(['Both groups exhibit a similar proportion of sessions with signficant '...
%             'phase-amplitude coupling, X(',num2str(df),')=',num2str(chi2stat),', p=',num2str(p)])
%     end
%     
%     HG_group1=Group1(Group1(:,5)==1,4);
%     HG_group2=Group2(Group2(:,5)==1,4);
%     
%     CDFplots(HG_group1,HG_group2,{'Sacc','PAE'},{'HG Modulation Index'},2)
%     
%     
%     if ~isempty(LG_group2) && ~isempty(LG_group1)
%         %test variance
%         h = vartest2(HG_group1,HG_group2);
%         
%         if h==0 %run t-test with equal variances assumed
%             [h,p,ci,HGstats] = ttest2(HG_group1,HG_group2,'Vartype','equal');
%         else %run t-test controlling for unequal variances
%             [h,p,ci,HGstats] = ttest2(HG_group1,HG_group2,'Vartype','unequal');
%         end
%         
%         HGmu_group1=nanmean(Group1(Group1(:,5)==1,4));
%         HGmu_group2=nanmean(Group2(Group2(:,5)==1,4));
%         
%         HGsem_group1=std(Group1(Group1(:,5)==1,4))/sqrt(size(Group1(Group1(:,5)==1,4),1));
%         HGsem_group2=std(Group2(Group2(:,5)==1,4))/sqrt(size(Group2(Group2(:,5)==1,4),1));
%         
%         disp(['High Gamma Magnitude of MI in area ',area{region},' session ',num2str(sess)])
%         if p<.05 && HGmu_group1<HGmu_group2
%             disp(['PAE rats exhibited higher mean MI, mean=',num2str(LGmu_group2),'sem=',num2str(LGsem_group2),' across the entire session'...
%                 'relative to Sacc controls, mean=',num2str(LGmu_group1),'sem=',num2str(LGsem_group1),'p=',num2str(p)])
%         elseif p<.05 && HGmu_group1>HGmu_group2
%             disp(['Sacc rats exhibited higher mean MI, mean=',num2str(HGmu_group1),'sem=',num2str(HGsem_group1),' across the entire session'...
%                 'relative to PAE rats, mean=',num2str(HGmu_group2),'sem=',num2str(HGsem_group2),'p=',num2str(p)])
%         else
%             disp(['There is no significant difference between MI between groups, p=',num2str(p)])
%         end
%     end
end
