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

%% RYAN DATA
% cd D:\Projects\PAE_PlaceCell\CFCLapResults

%% LAURA DATA 
 cd D:\Projects\PAE_PlaceCell\CFCLapResults

% rats=dir('*.mat');
% rats={rats.name};



for region=1:2
    %INITIALIZE DATA STRUCTURES
    
    tempResults.controls.sess{1}=[];
    tempResults.controls.sess{2}=[];
    tempResults.controls.sess{3}=[];
    tempResults.controls.sess{4}=[];
    
    tempResults.pae.sess{1}=[];
    tempResults.pae.sess{2}=[];
    tempResults.pae.sess{3}=[];
    tempResults.pae.sess{4}=[];
    
    tempResults.controls.laps.power.theta=[];
    tempResults.pae.laps.power.theta=[];
    tempResults.controls.laps.power.LG=[];
    tempResults.pae.laps.power.LG=[];
    tempResults.controls.laps.power.HG=[];
    tempResults.pae.laps.power.HG=[];
    tempResults.controls.laps.LGMI=[];
    tempResults.pae.laps.LGMI=[];
    tempResults.controls.laps.HGMI=[];
    tempResults.pae.laps.HGMI=[];
    tempResults.controls.laps.LGMI_thres=[];
    tempResults.pae.laps.LGMI_thres=[];
    tempResults.controls.laps.HGMI_thres=[];
    tempResults.pae.laps.HGMI_thres=[];
    
    rats=stoops{region};
    
    for i=1:length(rats)
        
        %Check if session exists
        if exist([rats{i},'.mat'],'file')
            cfc=load([rats{i},'.mat']);
        else
            continue
        end
        
        %Loop through sessions
        for sess=1:size(cfc.cfcLap,2)
            
            %For track sessions that also have the power field
            if contains(cfc.cfcLap{1,sess}.mazetype,["LinearTrack","linear_track","CircularTrack"]) && isfield(cfc.cfcLap{1,sess},'power')
                
                %Make sure there are at least 10 laps
                if size(cfc.cfcLap{1, 1}.power.rightLap.LG,2)<10 || size(cfc.cfcLap{1, 1}.power.leftLap.LG,2)<10
                    continue
                end
                
                %GET THETA AND GAMMA POWER FOR EACH LAP
                vars=fieldnames(cfc.cfcLap{1, sess}.power.rightLap.LG{1, 1}); %get varnames
                powVars=vars(contains(vars,'Pow')); %create index for power vars
                MIvars=vars(~contains(vars,{'Pow','_amps','_out'},'IgnoreCase',true));
                
                
                %Right Laps theta & low gamma
                for   lap=1:10 %iterate through first 10 laps
                    %power loop
                    for var=1:size(powVars,1) %iterate through vars
                        rightData(lap,var)=mean(getfield(cfc.cfcLap{1, 1}.power.rightLap.LG{1, lap},powVars{var}));
                    end
                    %MI loop
                    for var=1:size(MIvars,1) %iterate through vars
                        rightMI(lap,var)=getfield(cfc.cfcLap{1, 1}.power.rightLap.LG{1, lap},MIvars{var});
                    end
                end
                
                right_meanTheta=rightData(:,1);
                right_meanLG=rightData(:,2);
                
                rightLGMI_obtained=rightMI(:,1); %phase amplitude coupling measure - theta/low gamma
                rightLGMI_thres=rightMI(:,2);
                rightLGMI_sig=rightLGMI_obtained>rightLGMI_thres;
                
                %Right Laps high gamma
                for   lap=1:10 %iterate through first 10 laps
                    for var=1:size(powVars,1) %iterate through vars
                        rightData(lap,var)=mean(getfield(cfc.cfcLap{1, 1}.power.rightLap.HG{1, lap},powVars{var}));
                    end
                    %MI loop
                    for var=1:size(MIvars,1) %iterate through vars
                        rightMI(lap,var)=getfield(cfc.cfcLap{1, 1}.power.rightLap.HG{1, lap},MIvars{var});
                    end
                end
                
                right_meanHG=rightData(:,2);
                
                rightHGMI_obtained=rightMI(:,1); %phase amplitude coupling measure - theta/high gamma
                rightHGMI_thres=rightMI(:,2);
                rightHGMI_sig=rightHGMI_obtained>rightHGMI_thres;
                

                clear vars powVars lap 
                
                vars=fieldnames(cfc.cfcLap{1, sess}.power.leftLap.LG{1, 1}); %get varnames
                powVars=vars(contains(vars,'Pow')); %create index for power vars
                MIvars=vars(~contains(vars,{'Pow','_amps','_out'},'IgnoreCase',true));
                
                %Left Laps theta & low gamma
                for   lap=1:10%iterate through laps
                    % Power loop
                    for var=1:size(powVars,1) %iterate through vars
                        leftData(lap,var)=mean(getfield(cfc.cfcLap{1, 1}.power.leftLap.LG{1, lap},powVars{var}));
                    end
                    %MI loop
                    for var=1:size(MIvars,1) %iterate through vars
                        leftMI(lap,var)=getfield(cfc.cfcLap{1, 1}.power.leftLap.LG{1, lap},MIvars{var});
                    end
                    
                end
                
                left_meanTheta=leftData(:,1);
                left_meanLG=leftData(:,2);
       
                leftLGMI_obtained=leftMI(:,1); %phase amplitude coupling measure - theta/low gamma
                leftLGMI_thres=leftMI(:,2);
                leftLGMI_sig=leftLGMI_obtained>leftLGMI_thres;
                
                %Left Laps high gamma
                for   lap=1:10 %iterate through first 10 laps
                    for var=1:size(powVars,1) %iterate through vars
                        leftData(lap,var)=mean(getfield(cfc.cfcLap{1, 1}.power.leftLap.HG{1, lap},powVars{var}));
                    end
                    %MI loop
                    for var=1:size(MIvars,1) %iterate through vars
                        leftMI(lap,var)=getfield(cfc.cfcLap{1, 1}.power.leftLap.HG{1, lap},MIvars{var});
                    end
                end
                
                left_meanHG=leftData(:,2);
                
                leftHGMI_obtained=leftMI(:,1);
                leftHGMI_thres=leftMI(:,2);
                leftHGMI_sig=leftHGMI_obtained>leftHGMI_thres;
                
                %Build Data Structure for power & MI
                if contains(rats(i),control)
                    
                    tempResults.controls.laps.power.theta=[tempResults.controls.laps.power.theta;...
                        right_meanTheta'; left_meanTheta'];
                    
                    tempResults.controls.laps.power.LG=[tempResults.controls.laps.power.LG; ...
                        right_meanLG'; left_meanLG'];
                    
                    tempResults.controls.laps.power.HG=[tempResults.controls.laps.power.HG;...
                        right_meanHG'; left_meanHG'];
                    
                    tempResults.controls.laps.LGMI=[tempResults.controls.laps.LGMI; ...
                        rightLGMI_obtained'; leftLGMI_obtained'];
                    
                     tempResults.controls.laps.HGMI=[tempResults.controls.laps.LGMI; ...
                        rightHGMI_obtained'; leftHGMI_obtained'];
                    
                    tempResults.controls.laps.LGMI_thres=[tempResults.controls.laps.LGMI_thres; ...
                        rightLGMI_sig'; leftLGMI_sig'];
                    
                     tempResults.controls.laps.HGMI_thres=[tempResults.controls.laps.HGMI_thres; ...
                        rightHGMI_sig'; leftHGMI_sig'];
                    
                elseif contains(rats(i),pae)
                    
                     tempResults.pae.laps.power.theta=[tempResults.pae.laps.power.theta;...
                        right_meanTheta'; left_meanTheta'];
                    
                    tempResults.pae.laps.power.LG=[tempResults.pae.laps.power.LG; ...
                        right_meanLG'; left_meanLG'];
                    
                    tempResults.pae.laps.power.HG=[tempResults.pae.laps.power.HG;...
                        right_meanHG'; left_meanHG'];
                    
                     tempResults.pae.laps.LGMI=[tempResults.pae.laps.LGMI; ...
                        rightLGMI_obtained'; leftLGMI_obtained'];
                    
                     tempResults.pae.laps.HGMI=[tempResults.pae.laps.HGMI; ...
                        rightHGMI_obtained'; leftHGMI_obtained'];
                    
                    tempResults.pae.laps.LGMI_thres=[tempResults.pae.laps.LGMI_thres; ...
                        rightLGMI_sig'; leftLGMI_sig'];
                    
                     tempResults.pae.laps.HGMI_thres=[tempResults.pae.laps.HGMI_thres; ...
                        rightHGMI_sig'; leftHGMI_sig'];
                    
                end
                
                clear vars powVars lap leftHGMI_sig leftHGMI_thres leftHGMI_obtained...
                    quadHG_left lapHG_left leftLGMI_sig leftLGMI_thres leftLGMI_obtained...
                    quadLG_left lapLG_left...
                    quadtheta_left laptheta_left leftData leftMI
                
                
            else
              continue 
            end
            
        end
        clear cfc
    end
    
    
    clear MIvars var sess i
    
    %% RUN STATS

    %LAP ANALYSIS
    
    %power over first 10 laps
    
    disp('Lap analysis for theta power')
    
    Group1=tempResults.controls.laps.power.theta;
    Group2=tempResults.pae.laps.power.theta;
    
    disp(['Power over quadrants for theta'])
    
    %Combined laps
%     analyze_cfc(Group1,Group2,{'theta'});
    RL_anova(Group1,Group2,{'theta'})
    
    
    disp('Lap analysis for low gamma power')
    
    Group1=tempResults.controls.laps.power.LG;
    Group2=tempResults.pae.laps.power.LG;
    
    disp(['Power over quadrants for low gamma'])
    
    %Combined laps
    RL_anova(Group1,Group2,{'low gamma'});
    
     disp('Lap analysis for high gamma power')
    
    Group1=tempResults.controls.laps.power.HG;
    Group2=tempResults.pae.laps.power.HG;
    
    disp('Power over quadrants for High gamma')
    
    %Combined laps
    RL_anova(Group1,Group2,{'high gamma'});
    
    %MI OVER LAPS
    
    disp('Lap analysis low gamma amplitude theta phase coupling')
    
    Group1=tempResults.controls.laps.LGMI(:);
    Group2=tempResults.pae.laps.LGMI(:);
    
    g1Idx=logical(tempResults.controls.laps.LGMI_thres(:));
    g2Idx=logical(tempResults.pae.laps.LGMI_thres(:));
   
    disp('Proportion of Session with low gamma theta CFC')
    
    %low gamma
    group1_LGMI=sum(g1Idx)/size(Group1,1);
    group2_LGMI=sum(g2Idx)/size(Group2,1);
    
    %test proportion of sessions that exhibit phase-amplitude coupling
    [~,p, chi2stat,df] = prop_test([group1_LGMI group2_LGMI] , [1 1], 0)
    
    Group1=tempResults.controls.laps.HGMI(:);
    Group2=tempResults.pae.laps.HGMI(:);
    
    g1Idx=logical(tempResults.controls.laps.HGMI_thres(:));
    g2Idx=logical(tempResults.pae.laps.HGMI_thres(:));
   
    disp('Proportion of Session with low gamma theta CFC')
    
    %low gamma
    group1_HGMI=sum(g1Idx)/size(Group1,1);
    group2_HGMI=sum(g2Idx)/size(Group2,1);
    
    %test proportion of sessions that exhibit phase-amplitude coupling
    [~,p, chi2stat,df] = prop_test([group1_HGMI group2_HGMI] , [1 1], 0);
    
    %TRY TO LOOK AT MI OVER TRIALS
    Group1=tempResults.controls.laps.LGMI;
    Group2=tempResults.pae.laps.LGMI;
    
%     g1Idx=logical(tempResults.controls.laps.LGMI_thres);
%     g2Idx=logical(tempResults.pae.laps.LGMI_thres);
%     
%     Group1(~g1Idx)=NaN;
%     Group2(~g2Idx)=NaN;
%     
%     Group1=Group1(sum(isnan(Group1),2)==size(Group1,2),:);
%         Group2=Group2(sum(isnan(Group2),2)==size(Group2,2),:);

%     [p,tbl] = friedman(Group1)
    
    RL_anova(Group1,Group2,{'theta-LG MI'});

    
%     meanG1=nanmedian(Group1);
%     meanG2=nanmedian(Group2);
    
    
      %TRY TO LOOK AT MI OVER TRIALS
    Group1=tempResults.controls.laps.HGMI;
    Group2=tempResults.pae.laps.HGMI;
        RL_anova(Group1,Group2,{'theta-HG MI'});

    
    g1Idx=logical(tempResults.controls.laps.HGMI_thres);
    g2Idx=logical(tempResults.pae.laps.HGMI_thres);
    
    Group1(~g1Idx)=NaN;
    Group2(~g2Idx)=NaN;
   
end
