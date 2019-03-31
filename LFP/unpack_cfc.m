clear
clc
close all

control={'RH13','RH14','LS21','LS23','LE2821','LE2823','LEM3116','LEM3120'};
pae={'RH11','RH16','LS17','LS19','LE2813'};


%% SPLIT BY REGION
% load metadata files and extract region info
cd D:\Projects\PAE_PlaceCell\AnimalMetadata

rats=dir('*.mat');
rats={rats.name};
sess_region=[];
sessionid=[];
% mainpath='D:\Projects\PAE_PlaceCell\ProcessedData\';
mainpath=[];

for i=1:length(rats)
    load(rats{i})
    sess=fieldnames(AnimalMetadata.RecordingLogs);
    for s=1:length(sess)
        sess_region=[sess_region;{AnimalMetadata.AnimalName,sess{s},AnimalMetadata.RecordingLogs.(sess{s}).RecordingArea}];
        sessionid=[sessionid;{[mainpath,AnimalMetadata.AnimalName,'_',sess{s}]}];
    end
end

% create idx
ca1idx=strcmp(sess_region(:,3), 'ca1');
ca3idx=strcmp(sess_region(:,3), 'ca3');

ca1=sessionid(ca1idx);
ca3=sessionid(ca3idx);

stoops{1}=ca1;
stoops{2}=ca3;
%%
cd D:\Projects\PAE_PlaceCell\CFCResults

% rats=dir('*.mat');
% rats={rats.name};

%INITIALIZE DATA STRUCTURES

tempResults.controls.sess{1}=[];
tempResults.controls.sess{2}=[];
tempResults.controls.sess{3}=[];
tempResults.controls.sess{4}=[];

tempResults.pae.sess{1}=[];
tempResults.pae.sess{2}=[];
tempResults.pae.sess{3}=[];
tempResults.pae.sess{4}=[];

tempResults.controls.right.Quad=[];
tempResults.controls.right.Laps=[];
tempResults.controls.right.QuadMI=[];

tempResults.controls.left.Quad=[];
tempResults.controls.left.Laps=[];
tempResults.controls.left.QuadMI=[];

tempResults.pae.right.Quad=[];
tempResults.pae.right.Laps=[];
tempResults.pae.right.QuadMI=[];

tempResults.pae.left.Quad=[];
tempResults.pae.left.Laps=[];
tempResults.pae.left.QuadMI=[];

for region=1:2
    rats=stoops{region};
    for i=1:length(rats)
        if exist([rats{i},'.mat'],'file')
            cfc=load(rats{i});
        else
            continue
        end
        
        

        
        for sess=1:size(cfc.cfc,2)

            if contains(cfc.cfc{1,sess}.mazetype,["LinearTrack","linear_track","CircularTrack"]) && isfield(cfc.cfc{1,sess},'power')
                
                if contains(rats(i),control)
                    tempResults.controls.RatID{i,1}=rats{i};
                elseif contains(rats(i),pae)
                    tempResults.pae.RatID{i,1}=rats{i};
                end
                
                if size(cfc.cfc{1, 1}.power.rightLap.LG,2)<10 || size(cfc.cfc{1, 1}.power.leftLap.LG,2)<10
                    continue
                end
                
                %GET THETA AND GAMMA POWER FOR EACH LAP
                vars=fieldnames(cfc.cfc{1, sess}.power.rightLap.LG{1, 1}); %get varnames
                powVars=vars(contains(vars,'Pow')); %create index for power vars
                MIvars=vars(~contains(vars,{'Pow','_amps','_out'},'IgnoreCase',true));
                
                
                %Right Laps theta & low gamma
                for   lap=1:10 %iterate through first 10 laps
                    %power loop
                    for var=1:size(powVars,1) %iterate through vars
                        rightData(lap,var)=mean(getfield(cfc.cfc{1, 1}.power.rightLap.LG{1, lap},powVars{var}));
                    end
                    %MI loop
                    for var=1:size(MIvars,1) %iterate through vars
                        rightMI(lap,var)=getfield(cfc.cfc{1, 1}.power.rightLap.LG{1, lap},MIvars{var});
                    end
                end
                
                laptheta_right=mean(rightData(:,[1,3,5,7]),2); %theta collapse across lap
                quadtheta_right=mean(rightData(:,[1,3,5,7]));%theta collapse across quad
                
                lapLG_right=mean(rightData(:,[2,4,6,8]),2); %low gamma collapsed across lap
                quadLG_right=mean(rightData(:,[2,4,6,8])); %low gamma collapsed across quad
                
                rightLGMI_obtained=rightMI(:,[1,3,5,7]); %phase amplitude coupling measure - theta/low gamma
                rightLGMI_thres=rightMI(:,[2,4,6,8]);
                rightLGMI_sig=rightLGMI_obtained>rightLGMI_thres;
                
                %Right Laps high gamma
                for   lap=1:10 %iterate through first 10 laps
                    for var=1:size(powVars,1) %iterate through vars
                        rightData(lap,var)=mean(getfield(cfc.cfc{1, 1}.power.rightLap.HG{1, lap},powVars{var}));
                    end
                    %MI loop
                    for var=1:size(MIvars,1) %iterate through vars
                        rightMI(lap,var)=getfield(cfc.cfc{1, 1}.power.rightLap.HG{1, lap},MIvars{var});
                    end
                end
                
                quadHG_right=mean(rightData(:,[2,4,6,8])); %high gamma collapsed across lap
                lapHG_right=mean(rightData(:,[2,4,6,8]),2); %high gamma collapsed across quad
                
                rightHGMI_obtained=rightMI(:,[1,3,5,7]); %phase amplitude coupling measure - theta/high gamma
                rightHGMI_thres=rightMI(:,[2,4,6,8]);
                rightHGMI_sig=rightHGMI_obtained>rightHGMI_thres;
                
                
                %Build Data Structure for power
                if contains(rats(i),control)
                    
                    tempResults.controls.right.Quad=[tempResults.controls.right.Quad; quadtheta_right...
                        quadLG_right quadHG_right];
                    tempResults.controls.right.Laps=[tempResults.controls.right.Laps; laptheta_right'...
                        lapLG_right' lapHG_right'];
                    tempResults.controls.right.QuadMI=[tempResults.controls.right.QuadMI;...
                        rightLGMI_obtained rightLGMI_sig rightHGMI_obtained rightHGMI_sig];
                    
                elseif contains(rats(i),pae)
                    
                    tempResults.pae.right.Quad=[tempResults.pae.right.Quad; quadtheta_right...
                        quadLG_right quadHG_right];
                    tempResults.pae.right.Laps=[tempResults.pae.right.Laps; laptheta_right'...
                        lapLG_right' lapHG_right'];
                    tempResults.pae.right.QuadMI=[tempResults.pae.right.QuadMI;...
                        rightLGMI_obtained rightLGMI_sig rightHGMI_obtained rightHGMI_sig];
                    
                end
                
                clear vars powVars lap rightHGMI_sig rightHGMI_thres rightHGMI_obtained...
                    quadHG_right lapHG_right rightLGMI_sig rightLGMI_thres rightLGMI_obtained...
                    quadLG_right lapLG_right...
                    quadtheta_right laptheta_right rightData rightMI
                
                vars=fieldnames(cfc.cfc{1, sess}.power.leftLap.LG{1, 1}); %get varnames
                powVars=vars(contains(vars,'Pow')); %create index for power vars
                MIvars=vars(~contains(vars,{'Pow','_amps','_out'},'IgnoreCase',true));
                
                %Left Laps theta & low gamma
                for   lap=1:10%iterate through laps
                    % Power loop
                    for var=1:size(powVars,1) %iterate through vars
                        leftData(lap,var)=mean(getfield(cfc.cfc{1, 1}.power.leftLap.LG{1, lap},powVars{var}));
                    end
                    %MI loop
                    for var=1:size(MIvars,1) %iterate through vars
                        leftMI(lap,var)=getfield(cfc.cfc{1, 1}.power.leftLap.LG{1, lap},MIvars{var});
                    end
                    
                end
                
                quadtheta_left=mean(leftData(:,[1,3,5,7])); %theta collapsed across quad
                laptheta_left=mean(leftData(:,[1,3,5,7]),2);%theta collapsed across lap
                
                lapLG_left=mean(leftData(:,[2,4,6,8]),2); %low gamma collapsed across lap
                quadLG_left=mean(leftData(:,[2,4,6,8])); %low gamma collapsed across quad
                
                leftLGMI_obtained=leftMI(:,[1,3,5,7]); %phase amplitude coupling measure - theta/low gamma
                leftLGMI_thres=leftMI(:,[2,4,6,8]);
                leftLGMI_sig=leftLGMI_obtained>leftLGMI_thres;
                
                %Left Laps high gamma
                for   lap=1:10 %iterate through first 10 laps
                    for var=1:size(powVars,1) %iterate through vars
                        leftData(lap,var)=mean(getfield(cfc.cfc{1, 1}.power.leftLap.HG{1, lap},powVars{var}));
                    end
                    %MI loop
                    for var=1:size(MIvars,1) %iterate through vars
                        leftMI(lap,var)=getfield(cfc.cfc{1, 1}.power.leftLap.HG{1, lap},MIvars{var});
                    end
                end
                
                quadHG_left=mean(leftData(:,[2,4,6,8]));
                lapHG_left=mean(leftData(:,[2,4,6,8]),2);
                
                leftHGMI_obtained=leftMI(:,[1,3,5,7]);
                leftHGMI_thres=leftMI(:,[2,4,6,8]);
                leftHGMI_sig=leftHGMI_obtained>leftHGMI_thres;
                
                %Build Data Structure for power
                if contains(rats(i),control)
                    
                    
                    tempResults.controls.left.Quad=[tempResults.controls.left.Quad; quadtheta_left...
                        quadLG_left quadLG_left];
                    tempResults.controls.left.Laps=[tempResults.controls.left.Laps; laptheta_left'...
                        lapLG_left' lapLG_left'];
                    tempResults.controls.left.QuadMI=[tempResults.controls.left.QuadMI;...
                        leftLGMI_obtained leftLGMI_sig leftHGMI_obtained leftHGMI_sig];
                    
                elseif contains(rats(i),pae)
                    
                    tempResults.pae.left.Quad=[tempResults.pae.left.Quad; quadtheta_left...
                        quadLG_left quadHG_left];
                    tempResults.pae.left.Laps=[tempResults.pae.left.Laps; laptheta_left'...
                        lapLG_left' lapHG_left'];
                    tempResults.pae.left.QuadMI=[tempResults.pae.left.QuadMI;...
                        leftLGMI_obtained leftLGMI_sig leftHGMI_obtained leftHGMI_sig];
                    
                end
                
                clear vars powVars lap leftHGMI_sig leftHGMI_thres leftHGMI_obtained...
                    quadHG_left lapHG_left leftLGMI_sig leftLGMI_thres leftLGMI_obtained...
                    quadLG_left lapLG_left...
                    quadtheta_left laptheta_left leftData leftMI
                
                peakTheta=cfc.cfc{1, sess}.theta_ch.f_peak;
                LG_MI=cfc.cfc{1, sess}.LGtortCFC.modindex;
                LG_MI_thres=cfc.cfc{1, sess}.LGtortCFC.MIthres;
                LG_MI_sig=LG_MI>LG_MI_thres;
                
                
                HG_MI=cfc.cfc{1, sess}.HGtortCFC.modindex;
                HG_MI_thres=cfc.cfc{1, sess}.HGtortCFC.MIthres;
                HG_MI_sig=HG_MI>HG_MI_thres;
                
                %Build Data Structure for whole session MI
                if contains(rats(i),control)
                    
                    tempResults.controls.sess{sess}=[tempResults.controls.sess{sess}; ...
                        peakTheta LG_MI LG_MI_sig HG_MI HG_MI_sig];
                    
                elseif contains(rats(i),pae)
                    
                    tempResults.pae.sess{sess}=[tempResults.pae.sess{sess}; ...
                        peakTheta LG_MI LG_MI_sig HG_MI HG_MI_sig];
                    
                end
                
            else
                
                peakTheta=cfc.cfc{1, sess}.theta_ch.f_peak;
                LG_MI=cfc.cfc{1, sess}.LGtortCFC.modindex;
                LG_MI_thres=cfc.cfc{1, sess}.LGtortCFC.MIthres;
                LG_MI_sig=LG_MI>LG_MI_thres;
                
                
                HG_MI=cfc.cfc{1, sess}.HGtortCFC.modindex;
                HG_MI_thres=cfc.cfc{1, sess}.HGtortCFC.MIthres;
                HG_MI_sig=HG_MI>HG_MI_thres;
                
                %Build Data Structure for whole session MI
                if contains(rats(i),control)
                    
                    tempResults.controls.sess{sess}=[tempResults.controls.sess{sess}; ...
                        peakTheta LG_MI LG_MI_sig HG_MI HG_MI_sig];
                    
                elseif contains(rats(i),pae)
                    
                    tempResults.pae.sess{sess}=[tempResults.pae.sess{sess}; ...
                        peakTheta LG_MI LG_MI_sig HG_MI HG_MI_sig];
                    
                end
                
                clear peakTheta LG_MI LG_MI_thres LG_MI_sig HG_MI HG_MI_thres HG_MI_sig
            end
            
        end
        clear cfc
    end
    
    
    clear MIvars var sess i
    
    %% RUN STATS
    
    %QUAD ANALYSIS
    
    %power
    
    disp('Quadrant analysis for theta, low gamma and high gamma power')
    
    Group1_right=tempResults.controls.right.Quad;
    Group1_left=tempResults.controls.left.Quad;
    
    Group2_right=tempResults.pae.right.Quad;
    Group2_left=tempResults.pae.left.Quad;
    
    %freq power over quadrants
    freqVars={'theta','low gamma','high gamma'};
    frequency=1;
    
    for freq=1:4:12
        disp(['Power over quadrants for ',freqVars(frequency)])
        
        Group1_right_freq=Group1_right(:,freq:freq+3);
        Group1_left_freq=Group1_left(:,freq:freq+3);
        Group1_combined=[Group1_right_freq; fliplr(Group1_left_freq)];
        
        Group2_right_freq=Group2_right(:,freq:freq+3);
        Group2_left_freq=Group2_left(:,freq:freq+3);
        Group2_combined=[Group2_right_freq; fliplr(Group2_left_freq)];
        
        %     %Right test
        %     disp('Group differences - right laps - quadrants')
        %     analyze_cfc(Group1_right_freq,Group2_right_freq,freqVars(frequency));
        %
        %     %Left test
        %     disp('Group differences - left laps - quadrants')
        %     analyze_cfc(Group1_left_freq,Group2_left_freq,freqVars(frequency));
        
        %Combined laps
        disp('Group differences - combined - quadrants')
        analyze_cfc(Group1_combined,Group2_combined,freqVars(frequency));
        
        frequency=frequency+1;
        
    end
    %
    % %power
    %
    % disp('Quadrant analysis for phase-amplitude coupling (MI)')
    %
    % Group1=[tempResults.controls.right.QuadMI; fli ]
    % Group2=[tempResults.pae.right.QuadMI;  ]
    %
    %
    % %low gamma
    % LG_propSig_group1=sum(Group1(:,3))/size(Group1(:,2),1);
    % LG_propSig_group2=sum(Group2(:,3))/size(Group2(:,2),1);
    %
    % Group2_right=tempResults.pae.right.Quad;
    % Group2_left=tempResults.pae.left.Quad;
    %
    % %freq power over quadrants
    % freqVars={'theta','low gamma','high gamma'};
    % frequency=1;
    % for freq=1:4:12
    %     disp(['Power over quadrants for ',freqVars(frequency)])
    %
    %     Group1_right_freq=Group1_right(:,freq:freq+3);
    %     Group1_left_freq=Group1_left(:,freq:freq+3);
    %     Group1_combined=[Group1_right_freq; fliplr(Group1_left_freq)];
    %
    %     Group2_right_freq=Group2_right(:,freq:freq+3);
    %     Group2_left_freq=Group2_left(:,freq:freq+3);
    %     Group2_combined=[Group2_right_freq; fliplr(Group2_left_freq)];
    %
    %     %Right test
    %     disp('Group differences - right laps - quadrants')
    %     analyze_cfc(Group1_right_freq,Group2_right_freq,freqVars(frequency));
    %
    %     %Left test
    %     disp('Group differences - left laps - quadrants')
    %     analyze_cfc(Group1_left_freq,Group2_left_freq,freqVars(frequency));
    %
    %     %Combined laps
    %     disp('Group differences - combined - quadrants')
    %     analyze_cfc(Group1_combined,Group2_combined,freqVars(frequency));
    %
    %     frequency=frequency+1;
    %
    % end
    %LAP ANALYSIS
    
    %FREQUENCY ANALYSIS
    
    % Frequency
    %
    % ThetaFreq = CDFplots(Group1(:,1),Group2(:,1),GroupNames,VarNames(1),2);
    
    %% WHOLE SESSION MI ANALYSIS
    
    for sess=1:3
        disp(['Session ', num2str(sess), ' Analysis'])
        VarNames={'Peak Theta Frequency','MI- Theta/Low Gamma','MI- Theta/Low Gamma'};
        GroupNames={'Sacc','PAE'};
        
        %session 1 is linear track
        Group1=tempResults.controls.sess{1, sess};
        Group2=tempResults.pae.sess{1, sess};
        
        %low gamma
        LG_propSig_group1=sum(Group1(:,3))/size(Group1(:,2),1);
        LG_propSig_group2=sum(Group2(:,3))/size(Group2(:,2),1);
        
        %test proportion of sessions that exhibit phase-amplitude coupling
        [~,p, chi2stat,df] = prop_test([LG_propSig_group1 LG_propSig_group2] , [1 1], 0);
        
        disp('Results for theta-low gamma phase-amplitude coupling:')
        if p<.05
            disp(['The total number of sessions that exhibit significant'...
                'phase-amplitude coupling is different between groups, X(',num2str(df),')=',num2str(chi2stat),', p=',num2str(p)])
        else
            disp(['Both groups exhibit a similar proportion of sessions with signficant '...
                'phase-amplitude coupling, X(',num2str(df),')=',num2str(chi2stat),', p=',num2str(p)])
        end
        
        %test magnitude difference of MI between groups
        LG_group1=Group1(Group1(:,3)==1,2);
        LG_group2=Group2(Group2(:,3)==1,2);
        
        %test variance
        if ~isempty(LG_group2) && ~isempty(LG_group1)
            
            h = vartest2(LG_group1,LG_group2);
            
            CDFplots(LG_group1,LG_group2,{'Sacc','PAE'},{'LG Modulation Index'},2)
            
            
            
            if h==0 %run t-test with equal variances assumed
                [~,p,~,LGstats] = ttest2(LG_group1,LG_group2,'Vartype','equal');
            else %run t-test controlling for unequal variances
                [~,p,~,LGstats] = ttest2(LG_group1,LG_group2,'Vartype','unequal');
            end
            
            LGmu_group1=nanmean(Group1(Group1(:,3)==1,2));
            LGmu_group2=nanmean(Group2(Group2(:,3)==1,2));
            
            LGsem_group1=std(Group1(Group1(:,3)==1,2))/sqrt(size(Group1(Group1(:,3)==1,2),1));
            LGsem_group2=std(Group2(Group2(:,3)==1,2))/sqrt(size(Group2(Group2(:,3)==1,2),1));
            
            if p<.05 && LGmu_group1<LGmu_group2
                disp(['PAE rats exhibited higher mean MI, mean=',num2str(LGmu_group2),...
                    'sem=',num2str(LGsem_group2),' across the entire session'...
                    'relative to Sacc controls, mean=',num2str(LGmu_group1),'sem=',num2str(LGsem_group1),'p=',num2str(p)])
            elseif p<.05 && LGmu_group1>LGmu_group2
                disp(['Sacc rats exhibited higher mean MI, mean=',num2str(LGmu_group1),...
                    'sem=',num2str(LGsem_group1),' across the entire session'...
                    'relative to PAE rats, mean=',num2str(LGmu_group2),'sem=',num2str(LGsem_group2),'p=',num2str(p)])
            else
                disp(['There is no significant difference between MI between groups, p=',num2str(p)])
            end
            
        end
        
        
        %high gamma
        HG_propSig_group1=sum(Group1(:,5))/size(Group1(:,4),1);
        HG_propSig_group2=sum(Group2(:,5))/size(Group2(:,4),1);
        
        %test proportion of sessions that exhibit phase-amplitude coupling
        [~,p, chi2stat,df] = prop_test([HG_propSig_group1 HG_propSig_group2] , [1 1], 0);
        
        disp('Results for theta-high gamma phase-amplitude coupling:')
        if p<.05
            disp(['The total number of sessions that exhibit significant'...
                'phase-amplitude coupling is different between groups,X(',num2str(df),')=',num2str(chi2stat),', p=',num2str(p)])
        else
            disp(['Both groups exhibit a similar proportion of sessions with signficant '...
                'phase-amplitude coupling, X(',num2str(df),')=',num2str(chi2stat),', p=',num2str(p)])
        end
        
        HG_group1=Group1(Group1(:,5)==1,4);
        HG_group2=Group2(Group2(:,5)==1,4);
        
                    CDFplots(HG_group1,HG_group2,{'Sacc','PAE'},{'HG Modulation Index'},2)

        
        if ~isempty(LG_group2) && ~isempty(LG_group1)
            %test variance
            h = vartest2(HG_group1,HG_group2);
            
            if h==0 %run t-test with equal variances assumed
                [h,p,ci,HGstats] = ttest2(HG_group1,HG_group2,'Vartype','equal');
            else %run t-test controlling for unequal variances
                [h,p,ci,HGstats] = ttest2(HG_group1,HG_group2,'Vartype','unequal');
            end
            
            HGmu_group1=nanmean(Group1(Group1(:,5)==1,4));
            HGmu_group2=nanmean(Group2(Group2(:,5)==1,4));
            
            HGsem_group1=std(Group1(Group1(:,5)==1,4))/sqrt(size(Group1(Group1(:,5)==1,4),1));
            HGsem_group2=std(Group2(Group2(:,5)==1,4))/sqrt(size(Group2(Group2(:,5)==1,4),1));
            
            if p<.05 && HGmu_group1<HGmu_group2
                disp(['PAE rats exhibited higher mean MI, mean=',num2str(LGmu_group2),'sem=',num2str(LGsem_group2),' across the entire session'...
                    'relative to Sacc controls, mean=',num2str(LGmu_group1),'sem=',num2str(LGsem_group1),'p=',num2str(p)])
            elseif p<.05 && HGmu_group1>HGmu_group2
                disp(['Sacc rats exhibited higher mean MI, mean=',num2str(HGmu_group1),'sem=',num2str(HGsem_group1),' across the entire session'...
                    'relative to PAE rats, mean=',num2str(HGmu_group2),'sem=',num2str(HGsem_group2),'p=',num2str(p)])
            else
                disp(['There is no significant difference between MI between groups, p=',num2str(p)])
            end
        end
    end
end
