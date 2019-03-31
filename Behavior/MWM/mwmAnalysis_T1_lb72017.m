% WaterMazeAnalysis analyzes watermaze data
% Input: Compiled XY coordinates from FIJI manual track plugin
%      Coordinates should be stored in an excel file named
%      T#_tracking i.e. T1_tracking (Timepoint 1 tracking file)
%
% Created by BC 2014, adapted for FIJI tracking by LB 2016-2017
clc;
close all;
clear;
tic
addpath(genpath('/Users/lauraberkowitz/Google Drive/MATLAB/chronux'))
addpath(genpath('C:\Users\Ben Clark''s Lab\Google Drive\MATLAB\chronux'))
addpath(genpath('C:\Users\Ben Clark''s Lab\Documents\MATLAB\MEX'))
addpath(genpath('C:\Users\Ben Clark''s Lab\Google Drive\MWM_Data'))
%Import Data
for j=1:13
%     cd 'C:\Users\Ben Clark''s Lab\Google Drive\MATLAB\MWM';
    close all
    tab={'Day1','Day2','Day3','Day4','Day5','Day6','Day7','Day8','Day9','Day10','Day11','Day12A', 'Day12B'}; %indicate types of tasks
    trackData=readtable('/Users/RyanHarvey/GoogleDrive/MatlabDir/BClarkToolbox/Analysis/MWM/T1_Tracking.xlsx','Sheet',tab{j});
    
    if j==1
        if isempty(trackData);
            continue
        end
        disp('Processing Day1')
        %Pool Parameters for Reference Memory Day 1
        T1x_Origin=220.444; T1y_Origin=255.5;T1x_platOrigin=151;T1y_platOrigin=334; T1x_max=443;T1x_min=0;T1y_max=449.461;T1y_min=48; theta = deg2rad(9.102);T1Fs=10; T1dsFs=10;
    end
    if j==2
        if isempty(trackData);
            continue
        end
        disp('Processing Day2')
        %Pool Parameters for Reference Memory Day 1
        T1x_Origin=237; T1y_Origin=271.5;T1x_platOrigin=168;T1y_platOrigin=355; T1x_max=448;T1x_min=26;T1y_max=465;T1y_min=66; theta = deg2rad(9.102);T1Fs=10; T1dsFs=10;
    end
    if j==3
        if isempty(trackData);
            continue
        end
        disp('Processing Day3')
        %Pool Parameters for Reference Memory Day 1
        T1x_Origin=218; T1y_Origin=274.5;T1x_platOrigin=145;T1y_platOrigin=361; T1x_max=432;T1x_min=4;T1y_max=472;T1y_min=74; theta = deg2rad(9.102);T1Fs=20; T1dsFs=10;
    end
    if j==4
        if isempty(trackData);
            continue
        end
        disp('Processing Day4')
        %Pool Parameters for Reference Memory Day 1
        T1x_Origin=220; T1y_Origin=274;T1x_platOrigin=155;T1y_platOrigin=361; T1x_max=429;T1x_min=11;T1y_max=473;T1y_min=74; theta = deg2rad(9.102);T1Fs=20; T1dsFs=10;
    end
    if j==5
        if isempty(trackData);
            continue
        end
        disp('Processing Day5')
        %Pool Parameters for Reference Memory Day 1
        T1x_Origin=234; T1y_Origin=279.5;T1x_platOrigin=167;T1y_platOrigin=362; T1x_max=450;T1x_min=17;T1y_max=477;T1y_min=81; theta = deg2rad(9.102);T1Fs=20; T1dsFs=10;
    end
    if j==6 %Probe Trial - Reference Memory
        if isempty(trackData);
            continue
        end
        disp('Processing Day6')
        %Pool Parameters for Reference Memory Probe
        T1x_Origin=240.5; T1y_Origin=280.784;T1x_platOrigin=151.931;T1y_platOrigin=364.785; T1x_max=449;T1x_min=30;T1y_max=449;T1y_min=66; theta = deg2rad(9.102);T1Fs=20; T1dsFs=10;
        %Create figure
        FigTot=figure('Name','Time Point 1 Probe Trial','NumberTitle','off'); %change title to fit data description
        Figf30=figure('Name','Time Point 1 Probe Trial: First 30s','NumberTitle','off'); %change title to fit data description
        Figl30=figure('Name','Time Point 1 Probe Trial: Last 30s','NumberTitle','off');%change title to fit data description
        Figf3=figure('Name','Time Point 1 Probe Trial: Initial Trajectory','NumberTitle','off')
        %Initilize subgroup matrices for heat maps/scatter box plots
        smoothTgM=[];
        smoothTgF=[];
        smoothWTM=[];
        smoothWTF=[];
        smoothTg=[];
        smoothWT=[];
    end
    if  j==7
        if isempty(trackData);
            continue
        end
        disp('Processing Day7')
        %Pool Parameters for M2P Day 1
        T1x_Origin=221.5; T1y_Origin=278.5;T1x_platOrigin=227;T1y_platOrigin=408; T1x_max=434;T1x_min=8;T1y_max=478;T1y_min=79; theta = deg2rad(0);T1Fs=10; T1dsFs=10;
        T1x_platRef=151.931;T1y_platRef=364.785; T1x_platPrior=227;T1y_platPrior=408;
    end
    if j==8
        if isempty(trackData);
            continue
        end
        disp('Processing Day8')
        %Pool Parameters for M2P Day 2
        T1x_platRef=151.931;T1y_platRef=364.785; T1x_platPrior=227;T1y_platPrior=408;
        T1x_Origin=241.5; T1y_Origin=272;T1x_platOrigin=280;T1y_platOrigin=167; T1x_max=457;T1x_min=23;T1y_max=470;T1y_min=74; theta = deg2rad(0);T1Fs=10; T1dsFs=10;
    end
    if j==9
        if isempty(trackData);
            continue
        end
        disp('Processing Day9')
        %Pool Parameters for M2P Day 3
        T1x_platRef=151.931;T1y_platRef=364.785; T1x_platPrior=280;T1y_platPrior=167;
        T1x_Origin=233.5; T1y_Origin=286.7;T1x_platOrigin=336;T1y_platOrigin=282; T1x_max=450;T1x_min=19;T1y_max=479;T1y_min=89; theta = deg2rad(0);T1Fs=10; T1dsFs=10;
    end
    if j==10
        if isempty(trackData);
            continue
        end
        disp('Processing Day10')
        %Pool Parameters for M2P Day 4
        T1x_platRef=151.931;T1y_platRef=364.785; T1x_platPrior=336;T1y_platPrior=282;
        T1x_Origin=247.5; T1y_Origin=281;T1x_platOrigin=178;T1y_platOrigin=189; T1x_max=464;T1x_min=29;T1y_max=478;T1y_min=81; theta = deg2rad(0);T1Fs=10; T1dsFs=10;
    end
    if j==11
        if isempty(trackData);
            continue
        end
        disp('Processing Day11')
        %Pool Parameters for M2P Day 5
        T1x_platRef=151.931;T1y_platRef=364.785; T1x_platPrior=178;T1y_platPrior=189;
        T1x_Origin=226.5; T1y_Origin=276;T1x_platOrigin=214;T1y_platOrigin=287; T1x_max=438;T1x_min=13;T1y_max=472;T1y_min=78; theta = deg2rad(0);T1Fs=10; T1dsFs=10;
    end
    if j==12
        if isempty(trackData);
            continue
        end
        disp('Processing Day12A')
        %Pool Parameters for Cued Group A
        T1x_Origin=240.5; T1y_Origin=280.784;T1x_platOrigin=151.931;T1y_platOrigin=364.785; T1x_max=449;T1x_min=30;T1y_max=449;T1y_min=66; theta = deg2rad(9.102);T1Fs=10; T1dsFs=10;
    end
     if j==13
        if isempty(trackData);
            continue
        end
        disp('Processing Day12B')
        %Pool Parameters for Cued Group B
        T1x_Origin=240.5; T1y_Origin=280.784;T1x_platOrigin=151.931;T1y_platOrigin=364.785; T1x_max=449;T1x_min=30;T1y_max=449;T1y_min=66; theta = deg2rad(9.102);T1Fs=10; T1dsFs=10;
    end
    

%%##### Convert Coordinates from pixels to cm ############################

%Pixel to cm converstion factor - 150 indicates pool diameter -
%____T1__________
T1x_conv=150/(T1x_max-T1x_min); %cm/pixel
T1y_conv=150/(T1y_max-T1y_min);

%Pool Boundaries in cm
%____T1__________
T1x_max_cm=(T1x_max*T1x_conv);
T1x_min_cm=(T1x_min*T1x_conv);
T1y_max_cm=(T1y_max*T1y_conv);
T1y_min_cm=(T1y_min*T1y_conv);

%Pool Origin in cm
T1x_Origin_cm=T1x_Origin*T1x_conv;
T1y_Origin_cm=T1y_Origin*T1y_conv;

%Platform Origin in cm
T1x_platOrigin_cm=T1x_platOrigin*T1x_conv;
T1y_platOrigin_cm=T1y_platOrigin*T1y_conv;

%Platform Origin transposed relative to pool origin of (0,0)
T1x_platOriginTranspose=T1x_platOrigin_cm-T1x_Origin_cm;
T1y_platOriginTranspose=T1y_platOrigin_cm-T1y_Origin_cm;

%Create squre platform (make lower left corner position
xpos=T1x_platOriginTranspose-8; %for a 16x16cm square platform
ypos=T1y_platOriginTranspose-8;

%Convert tracking data to cm
trackData{:,3}= (((trackData{:,3}))-(T1x_Origin))*T1x_conv;
trackData{:,4}=((trackData{:,4}-(T1y_Origin))*T1y_conv);

%=-=-=-=-=-===-=-=-=--== Rotate data points=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
plat=[T1x_platOriginTranspose;T1y_platOriginTranspose];
% create a matrix which will be used later in calculations with center of
% rotation at 0,0
centerPlat =[0; 0];
centerData = repmat([0; 0], 1, length(trackData{:,3}));
% define a counter-clockwise rotation matrix
% theta = deg2rad(9.102);     % dervied offline x=atan(o/a)
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
% do the rotation...
rotPlat = R*(plat - centerPlat) + centerPlat;% shift points in the plane so that the center of rotation is at the origin% apply the rotation about the origin
rotData = R*(trackData{:,3:4}' - centerData) + centerData;
% shift again so the origin goes back to the desired center of rotation
% pick out the vectors of rotated x- and y-data
x_rotated = rotData(1,:);
y_rotated = rotData(2,:);
trackData{:,3}=x_rotated';
trackData{:,4}=y_rotated';
T1x_platOriginRot= rotPlat(1,:);
T1y_platOriginRot= rotPlat(2,:);
%-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-==-Create Platform=-=-=-=-=-=-=-=-=-=-=-=-
%Create squre platform pos for rotated platform (make lower left corner position
xposRot=T1x_platOriginRot-8; %for a 16x16cm square platform
yposRot=T1y_platOriginRot-8;

%% Initialize variables to loop through xy coordinate data
% nameSplit=strsplit((name),'_'); %Identify timepoint via filename
[C, ia]=unique(trackData(:,1),'stable');
newC=reshape(C{:,1},[],1);
newC((sum(ismember(char(newC),'sbj'),2)~=3))=[];
%Pool Parameters
x=0;
y=0;
r=150/2; %150 is pool diameter
ang=0:0.01:2*pi;
xp=r*cos(ang);
yp=r*sin(ang);

%Plot Smoothed Data~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%Loop through each trial, process/plot data and compute output variables.
for i=1:size(newC,1) % %loop through #of participants
        try
            trialDsnuts=trackData(ia(i):(ia((i+1))-1),:);%Use unique index to process one session at a time.
        catch
        end
        if i==length(ia);
            trialDsnuts=trackData(ia(i):end,:);
        end
        trial=(strsplit(newC{i,:},'_'));
        trialData=reshape(trialDsnuts{:,2:4},[],3);%convert to double for processing with resample/runline
        trialData(isnan(trialData))=[];
        [p,q]=rat(T1dsFs/T1Fs);
        dsData=resample(trialData,p,q);
        
        if j== 1 || 2 || 3 || 4 ||5 || 6 && length(dsData)>600;%Probe session - shrink to 60sec duration
            dsData((length(dsData)-(((length(dsData)/10)-60)*10))+1:end,:)=[];
        end
        %Smooth Data
        x_smooth=runline(dsData(2:end,2),3,1); %Smooth data
        y_smooth=runline(dsData(2:end,3),3,1); %Smooth data
        
        %save trial data into .csv files by RH 07/2017
        Samples=1:size(x_smooth,1); TS=0:1/10:size(x_smooth,1);    TS=TS(1:size(x_smooth,1));
        id=trial;
        subject=strsplit(id{3},'j'); subject=subject(2);
        trialInfo={'Time_Point',num2str(id{1}),NaN,NaN;'Sample',num2str(size(x_smooth,1)),NaN,NaN;'ID',num2str(subject{1}),NaN,NaN;'trial',id{2},NaN,NaN};
        
        array_export=num2cell([Samples',TS',x_smooth,y_smooth]);
        
        exportsheet=[trialInfo;{'Samples','TS','X','Y'};array_export];
        
        filename=[pwd filesep 'Time_Point' num2str(id{1}) filesep strjoin({'day',num2str(j),trial{2},'_00',subject{1}},'') '.csv'];
        
        if i==1; mkdir([pwd filesep 'Time_Point' num2str(id{1})]); end
        
        fid=fopen(filename,'w');
        fprintf(fid,'%s, %s, %s, %s\n',exportsheet{1,:});
        fprintf(fid,'%s, %s, %s, %s\n',exportsheet{2,:});
        fprintf(fid,'%s, %s, %s, %s\n',exportsheet{3,:});
        fprintf(fid,'%s, %s, %s, %s\n',exportsheet{4,:});
        fprintf(fid,'%s, %s, %s, %s\n',exportsheet{5,:});
        fclose(fid);
        dlmwrite(filename,exportsheet(6:end,:),'-append');
        
        %Identify subgroups
        sbj=regexp(trial{3},'\d*', 'match');
        TgM={'1' '2' '3' '4' '5' '6' '7' '15'}; %subject IDs
        TgF={'8' '9' '10' '11' '12' '13' '14' '16'};
        WTM={'17' '18' '19' '20' '21' '22'};
        WTF={'23' '24' '25' '26' '27' '28'};
        Tg={'1' '2' '3' '4' '5' '6' '7' '15' '8' '9' '10' '11' '12' '13' '14' '16'};
        WT={'17' '18' '19' '20' '21' '22' '23' '24' '25' '26' '27' '28'};
        %-=-=-=-=-=-=-=-=-=-=-=-=-=-=-Obtain Outcome Measures Training-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        %Directional heading, path length, swim speed, circular statistics
        %for trials 9 - 20 (training)
        %Directional heading
        if  j==1 || j==2 || j==3 || j==4 || j==5; %Reference Memory
            %Plot data, add to figure whole session
            %Calculate cumulative proximity
            repTest=repmat([T1x_platOriginRot T1y_platOriginRot],([size(x_smooth,1) 1])); %initialize matrix for platform location
            proximity=[x_smooth y_smooth]- repTest(:,:); %calculate difference between path and platform location
            D(i,1)= (sum(sqrt((proximity(:,1).^2)+(proximity(:,2).^2))))-(sqrt(((x_smooth(1,1)-(T1x_platOriginRot))^2)+...
                (y_smooth(1,1)-(T1y_platOriginRot))^2)); %cumulative Proximity sans distance from start to plat origin.
            
            %Calculate Path length & Swim Speed
            sqrXDiff=(diff(x_smooth)).^2;
            sqrYDiff=(diff(y_smooth)).^2;
            pathDist=sqrt(sqrXDiff + sqrYDiff);
            PathLength(i,1)=sum(pathDist); %Total Path length
            swimSpeed(i,1)=(PathLength(i,1)/(size(x_smooth,1)/10)); %total swim speed
            
            %Calculate preference score total
            timeC=(sum(x_smooth(:)<=0 & y_smooth(:)>=0))/10; %calculate data points in target quadrant in seconds
            timeT=(sum(x_smooth(:)>=0 & y_smooth(:)>=0))/10;
            timeA=(sum(x_smooth(:)>=0 & y_smooth(:)<=0))/10;
            timeB=(sum(x_smooth(:)<=0 & y_smooth(:)<=0))/10;
            prefScoreTotal(i,1)=((timeT-timeA)+(timeT-timeB)+(timeT-timeC)/3);
            
            % Locate xy up to 45cm
            repTraject=repmat([x_smooth(1,1) y_smooth(1,1)],([(size(x_smooth,1)-1) 1]));
            trajectoryIndex=[x_smooth(2:end,1) y_smooth(2:end,1)]- repTraject(:,:);
            trajectoryDistance=sqrt((trajectoryIndex(:,1).^2)+(trajectoryIndex(:,2).^2));
            for idist=1:length(trajectoryDistance)
                if trajectoryDistance(idist)>=44.5; break; end;
            end
            hypotenuse=sqrt(((T1x_platOriginRot-x_smooth(1,1)).^2)+((T1y_platOriginRot-y_smooth(1,1)).^2));
            adjacent=sqrt(((T1x_platOriginRot-x_smooth(idist+1,1)).^2)+((T1y_platOriginRot-y_smooth(idist+1,1)).^2));
            trajectoryHeading(i,1)=acosd(adjacent/hypotenuse);
            trajectoryRadians(i,1)=acos(adjacent/hypotenuse);
            
            %Plot data, add to figure whole session
            %             set(0,'CurrentFigure',FigGroupA1);subplot(4,7,i)
            figure(i)
            plot(x+xp,y+yp);
            %     axis([-150 150]);
            title(newC(i),'FontSize', 12, 'Color','k','FontName','Helvetica');
            hold on
            rectangle('Position',[xposRot yposRot 16 16]);
            hold on
            plot(x_smooth,y_smooth);
            %plot quadrant boundaries
            n=4;
            tet=linspace(-pi,pi,n+1);
            xi=r*cos(tet)+x;
            yi=r*sin(tet)+y;
            for k=1:numel(xi)
                plot([x xi(k)],[y yi(k)])
                hold on
            end
            hold off
            
            %save paths to individual file
            saveas(figure(i),['C:\Users\Ben Clark''s Lab\Google Drive\MWM_Data\Analysis\AnalyzedData\Training\T1' filesep trial{3} filesep,newC{i},'.jpeg']);
            close all
            cd 'C:\Users\Ben Clark''s Lab\Google Drive\MATLAB\MWM'
            
        elseif j==6;
            %Plot data, add to figure whole session
            set(0,'CurrentFigure',FigTot);subplot(4,7,i)
            plot(x+xp,y+yp);
            %     axis([-150 150]);
            title(newC(i),'FontSize', 12, 'Color','k','FontName','Helvetica');
            hold on
            rectangle('Position',[xposRot yposRot 16 16]);
            hold on
            plot(x_smooth,y_smooth);
            %plot quadrant boundaries
            n=4;
            tet=linspace(-pi,pi,n+1);
            xi=r*cos(tet)+x;
            yi=r*sin(tet)+y;
            for k=1:numel(xi)
                plot([x xi(k)],[y yi(k)])
                hold on
            end
            hold off
            %Plot data, add to figure first 30 seconds
            set(0,'CurrentFigure',Figf30);subplot(4,7,i)
            plot(x+xp,y+yp);
            title(newC(i),'FontSize', 12, 'Color','k','FontName','Helvetica');
            hold on
            rectangle('Position',[xposRot yposRot 16 16]);
            hold on
            plot(x_smooth(1:300),y_smooth(1:300));
            %plot quadrant boundaries
            n=4;
            tet=linspace(-pi,pi,n+1);
            xi=r*cos(tet)+x;
            yi=r*sin(tet)+y;
            for k=1:numel(xi)
                plot([x xi(k)],[y yi(k)])
                hold on
            end
            hold off
            %Plot data, add to figure last 30 seconds
            set(0,'CurrentFigure',Figl30);subplot(4,7,i)
            plot(x+xp,y+yp);
            title(newC(i),'FontSize', 12, 'Color','k','FontName','Helvetica');
            hold on
            rectangle('Position',[xposRot yposRot 16 16]);
            hold on
            plot(x_smooth(301:end),y_smooth(301:end));
            %plot quadrant boundaries
            n=4;
            tet=linspace(-pi,pi,n+1);
            xi=r*cos(tet)+x;
            yi=r*sin(tet)+y;
            for k=1:numel(xi)
                plot([x xi(k)],[y yi(k)])
                hold on
            end
            hold off
            %Plot data, add to figure first 10sec
            set(0,'CurrentFigure',Figf3);subplot(4,7,i)
            plot(x+xp,y+yp);
            title(newC(i),'FontSize', 12, 'Color','k','FontName','Helvetica');
            hold on
            rectangle('Position',[xposRot yposRot 16 16]);
            hold on
            plot(x_smooth(1:100),y_smooth(1:100));
            %plot quadrant boundaries
            n=4;
            tet=linspace(-pi,pi,n+1);
            xi=r*cos(tet)+x;
            yi=r*sin(tet)+y;
            for k=1:numel(xi)
                plot([x xi(k)],[y yi(k)])
                hold on
            end
            hold off
            %extract subject data to indicate subgroup membership
            if ismember(sbj,TgM);
                %         catSubTgM = repmat(sbj,length(x_smooth),1);
                iTgM= [x_smooth y_smooth];
                smoothTgM=[smoothTgM; iTgM];
            elseif ismember(sbj,TgF);
                %         catSubTgF = repmat(sbj,length(x_smooth),1);
                iTgF=[x_smooth y_smooth];
                smoothTgF=[smoothTgF; iTgF];
            elseif ismember (sbj,WTM);
                %         catSubWTM =repmat(sbj,length(x_smooth),1);
                iWTM=[x_smooth y_smooth];
                smoothWTM=[smoothWTM; iWTM];
            else ismember(sbj,WTF);
                %         catSubWTF =repmat(sbj,length(x_smooth),1);
                iWTF=[x_smooth y_smooth];
                smoothWTF=[smoothWTF; iWTF];
            end
            
            %Compile data for heat maps between groups for probe data only - first 30s i.e. first 300
            %rows.
            if ismember(sbj,Tg) && strcmp('T11',trial{2});
                %         catSubTgM = repmat(sbj,length(x_smooth),1);
                iTg= [x_smooth(1:300,:) y_smooth(1:300,:)]; %concatenate first 30s
                smoothTg=[smoothTg; iTg];
            elseif ismember(sbj,WT) && strcmp('T11',trial{2});
                %         catSubTgF = repmat(sbj,length(x_smooth),1);
                iWT=[x_smooth(1:300,:) y_smooth(1:300,:)];%concatenate first 30s
                smoothWT=[smoothWT; iWT];
            end
            
            %~~~~~~~~~~~~~~~~~~~  Calculate Outcome Measures  ~~~~~~~~~~~~~~~~~~~~~~~
            %Calculate average proximity
            repTest=repmat([T1x_platOriginRot T1y_platOriginRot],([size(x_smooth,1) 1])); %initialize matrix for platform location
            proximity=[x_smooth y_smooth]- repTest(:,:); %calculate difference between path and platform location
            ProbeD(i,1)= mean(sqrt((proximity(:,1).^2)+(proximity(:,2).^2))); %Average Proximity
            ProbeDf30(i,1)=mean(sqrt((proximity(1:300,1).^2)+(proximity(1:300,2).^2)));
            ProbeDl30(i,1)=mean(sqrt((proximity(301:end,1).^2)+(proximity(301:end,2).^2)));
            
            %Create path differences between points
            sqrXDiff=(diff(x_smooth)).^2;
            sqrYDiff=(diff(y_smooth)).^2;
            pathDist=sqrt(sqrXDiff + sqrYDiff);
            ProbepathLength(i,1)=sum(pathDist); %Total Path length
            ProbeswimSpeed(i,1)=(ProbepathLength(i,1)/(size(x_smooth,1)/10));
            
            %Calculate initial trajectory
            % Locate xy up to 45cm
            repTraject=repmat([x_smooth(1,1) y_smooth(1,1)],([(size(x_smooth,1)-1) 1]));
            trajectoryIndex=[x_smooth(2:end,1) y_smooth(2:end,1)]- repTraject(:,:);
            trajectoryDistance=sqrt((trajectoryIndex(:,1).^2)+(trajectoryIndex(:,2).^2));
            for idist=1:length(trajectoryDistance)
                if trajectoryDistance(idist)>=44.5; break; end;
            end
            hypotenuse=sqrt(((T1x_platOriginRot-x_smooth(1,1)).^2)+((T1y_platOriginRot-y_smooth(1,1)).^2));
            adjacent=sqrt(((T1x_platOriginRot-x_smooth(idist+1,1)).^2)+((T1y_platOriginRot-y_smooth(idist+1,1)).^2));
            ProbetrajectoryHeading(i,1)=acosd(adjacent/hypotenuse);
            ProbetrajectoryRadians(i,1)=acos(adjacent/hypotenuse);
            
            %Calculate time in platform area - Quadrant T
            timePlatT(i,1)=(sum(x_smooth>=(T1x_platOriginRot-8) & x_smooth<=(T1x_platOriginRot+8) &...
                y_smooth>=(T1y_platOriginRot-8) & y_smooth<=(T1y_platOriginRot+8)))/10;
            
            %Calculate time in platform area - Quadrant A
            timePlatA(i,1)=(sum(x_smooth>=(-1.*(T1x_platOriginRot))-8 & x_smooth<=(-1.*(T1x_platOriginRot))+8 &...
                y_smooth>=T1y_platOriginRot-8 & y_smooth<=(T1y_platOriginRot+8)))/10;
            
            %Calculate time in platform area - Quadrant B
            timePlatB(i,1)=(sum(x_smooth>=(-1.*(T1x_platOriginRot))-8 & x_smooth<=(-1.*(T1x_platOriginRot))+8 &...
                y_smooth>=(-1.*(T1y_platOriginRot))-8 & y_smooth<=(-1.*(T1y_platOriginRot))+8))/10;
            
            %Calculate time in platform area - Quadrant C
            timePlatC(i,1)=(sum(x_smooth>=T1x_platOriginRot-8 & x_smooth<=T1x_platOriginRot+8 &...
                y_smooth>=(-1.*(T1y_platOriginRot))-8 & y_smooth<=(-1.*(T1y_platOriginRot))+8))/10;
            
            %Calculate preference score first 30s
            timeC=(sum(x_smooth(1:300)<=0 & y_smooth(1:300)>=0))/10; %calculate data points in target quadrant in seconds
            timeT=(sum(x_smooth(1:300)>=0 & y_smooth(1:300)>=0))/10;
            timeA=(sum(x_smooth(1:300)>=0 & y_smooth(1:300)<=0))/10;
            timeB=(sum(x_smooth(1:300)<=0 & y_smooth(1:300)<=0))/10;
            opQuadFirst(i,1)=timeB;
            prefScoreFirst(i,1)=((timeT-timeA)+(timeT-timeB)+(timeT-timeC)/3);
            
            %Calculate preference score last 30s
            timeC=(sum(x_smooth(301:end)<=0 & y_smooth(301:end)>=0))/10; %calculate data points in target quadrant in seconds
            timeT=(sum(x_smooth(301:end)>=0 & y_smooth(301:end)>=0))/10;
            timeA=(sum(x_smooth(301:end)>=0 & y_smooth(301:end)<=0))/10;
            timeB=(sum(x_smooth(301:end)<=0 & y_smooth(301:end)<=0))/10;
            opQuadLast(i,1)=timeB;
            prefScoreLast(i,1)=((timeT-timeA)+(timeT-timeB)+(timeT-timeC)/3);
            
            %Calculate preference score total
            timeC=(sum(x_smooth(:)<=0 & y_smooth(:)>=0))/10; %calculate data points in target quadrant in seconds
            timeT=(sum(x_smooth(:)>=0 & y_smooth(:)>=0))/10;
            timeA=(sum(x_smooth(:)>=0 & y_smooth(:)<=0))/10;
            timeB=(sum(x_smooth(:)<=0 & y_smooth(:)<=0))/10;
            opQuadTotal(i,1)=timeB;
            prefScoreTotal(i,1)=((timeT-timeA)+(timeT-timeB)+(timeT-timeC)/3);
            
        elseif strcmp('T11',trial{2})==0 && j==7 || j==8 || j==9 || j==10 || j==11; %M2P
            %Plot data, add to figure whole session
            %Calculate cumulative proximity
            repTest=repmat([T1x_platOriginRot T1y_platOriginRot],([size(x_smooth,1) 1])); %initialize matrix for platform location
            proximity=[x_smooth y_smooth]- repTest(:,:); %calculate difference between path and platform location
            D(i,1)= (sum(sqrt((proximity(:,1).^2)+(proximity(:,2).^2))))-(sqrt(((x_smooth(1,1)-(T1x_platOriginRot))^2)+...
                (y_smooth(1,1)-(T1y_platOriginRot))^2)); %cumulative Proximity sans distance from start to plat origin.
            
            %Calculate cumulative proximity to prior platform location
            repTest=repmat([((T1x_platPrior*T1x_conv)-T1x_Origin_cm) ((T1y_platPrior*T1y_conv)-T1y_Origin_cm)],([size(x_smooth,1) 1])); %initialize matrix for platform location
            proximity=[x_smooth y_smooth]- repTest(:,:); %calculate difference between path and platform location
            Dprior(i,1)= (sum(sqrt((proximity(:,1).^2)+(proximity(:,2).^2))))-(sqrt(((x_smooth(1,1)-((T1x_platPrior*T1x_conv)-T1x_Origin_cm))^2)+...
                (y_smooth(1,1)-((T1y_platPrior*T1y_conv)-T1y_Origin_cm))^2)); %cumulative Proximity sans distance from start to plat origin.
            
            %Calculate cumulative proximity to refernence memory platform
            %location
            repTest=repmat([((T1x_platRef*T1x_conv)-T1x_Origin_cm) ((T1y_platRef*T1y_conv)-T1y_Origin_cm)],([size(x_smooth,1) 1])); %initialize matrix for platform location
            proximity=[x_smooth y_smooth]- repTest(:,:); %calculate difference between path and platform location
            Dreference(i,1)= (sum(sqrt((proximity(:,1).^2)+(proximity(:,2).^2))))-(sqrt(((x_smooth(1,1)-((T1x_platRef*T1x_conv)-T1x_Origin_cm))^2)+...
                (y_smooth(1,1)-((T1y_platRef*T1y_conv)-T1y_Origin_cm))^2)); %cumulative Proximity sans distance from start to plat origin.
            
            %Calculate Path length & Swim Speed
            sqrXDiff=(diff(x_smooth)).^2;
            sqrYDiff=(diff(y_smooth)).^2;
            pathDist=sqrt(sqrXDiff + sqrYDiff);
            PathLength(i,1)=sum(pathDist); %Total Path length
            swimSpeed(i,1)=(PathLength(i,1)/(size(x_smooth,1)/10)); %total swim speed
            
            % Locate xy up to 45cm
            repTraject=repmat([x_smooth(1,1) y_smooth(1,1)],([(size(x_smooth,1)-1) 1]));
            trajectoryIndex=[x_smooth(2:end,1) y_smooth(2:end,1)]- repTraject(:,:);
            trajectoryDistance=sqrt((trajectoryIndex(:,1).^2)+(trajectoryIndex(:,2).^2));
            for idist=1:length(trajectoryDistance)
                if trajectoryDistance(idist)>=44.5; break; end;
            end
            hypotenuse=sqrt(((T1x_platOriginRot-x_smooth(1,1)).^2)+((T1y_platOriginRot-y_smooth(1,1)).^2));
            adjacent=sqrt(((T1x_platOriginRot-x_smooth(idist+1,1)).^2)+((T1y_platOriginRot-y_smooth(idist+1,1)).^2));
            trajectoryHeading(i,1)=acosd(adjacent/hypotenuse);
            trajectoryRadians(i,1)=acos(adjacent/hypotenuse);
            
            %Plot data, add to figure whole session
            %             set(0,'CurrentFigure',FigGroupA1);subplot(4,7,i)
            figure(i)
            plot(x+xp,y+yp);
            %     axis([-150 150]);
            title(newC(i),'FontSize', 12, 'Color','k','FontName','Helvetica');
            hold on
            rectangle('Position',[xposRot yposRot 16 16]);
            hold on
            plot(x_smooth,y_smooth);
            hold off
            %save paths to individual file
            saveas(figure(i),['C:\Users\Ben Clark''s Lab\Google Drive\MWM_Data\Analysis\AnalyzedData\M2P\T1' filesep trial{3} filesep,newC{i},'.jpeg']);
            close all
            cd 'C:\Users\Ben Clark''s Lab\Google Drive\MATLAB\MWM'
            
        elseif j==12 || j==13;
            %Plot data, add to figure whole session
            figure(i)
            plot(x+xp,y+yp);
            %     axis([-150 150]);
            title(newC(i),'FontSize', 12, 'Color','k','FontName','Helvetica');
            hold on
            rectangle('Position',[xposRot yposRot 16 16]);
            hold on
            plot(x_smooth,y_smooth);
            %plot quadrant boundaries
            n=4;
            tet=linspace(-pi,pi,n+1);
            xi=r*cos(tet)+x;
            yi=r*sin(tet)+y;
            for k=1:numel(xi)
                plot([x xi(k)],[y yi(k)])
                hold on
            end
            hold off
            saveas(figure(i),['C:\Users\Ben Clark''s Lab\Google Drive\MWM_Data\Analysis\AnalyzedData\Cued\T1' filesep trial{3} filesep,newC{i},'.jpeg']);
            close all
            cd 'C:\Users\Ben Clark''s Lab\Google Drive\MATLAB\MWM'
            %Calculate cumulative proximity
            repTest=repmat([T1x_platOriginRot T1y_platOriginRot],([size(x_smooth,1) 1])); %initialize matrix for platform location
            proximity=[x_smooth y_smooth]- repTest(:,:); %calculate difference between path and platform location
            D(i,1)= (sum(sqrt((proximity(:,1).^2)+(proximity(:,2).^2))))-(sqrt(((x_smooth(1,1)-(T1x_platOriginRot))^2)+...
                (y_smooth(1,1)-(T1y_platOriginRot))^2)); %cumulative Proximity sans distance from start to plat origin.
            
            %Calculate Path length & Swim Speed
            sqrXDiff=(diff(x_smooth)).^2;
            sqrYDiff=(diff(y_smooth)).^2;
            pathDist=sqrt(sqrXDiff + sqrYDiff);
            PathLength(i,1)=sum(pathDist); %Total Path length
            swimSpeed(i,1)=(PathLength(i,1)/(size(x_smooth,1)/10)); %total swim speed
            
            % Locate xy up to 45cm
            repTraject=repmat([x_smooth(1,1) y_smooth(1,1)],([(size(x_smooth,1)-1) 1]));
            trajectoryIndex=[x_smooth(2:end,1) y_smooth(2:end,1)]- repTraject(:,:);
            trajectoryDistance=sqrt((trajectoryIndex(:,1).^2)+(trajectoryIndex(:,2).^2));
            for idist=1:length(trajectoryDistance)
                if trajectoryDistance(idist)>=44.5; break; end;
            end
            hypotenuse=sqrt(((T1x_platOriginRot-x_smooth(1,1)).^2)+((T1y_platOriginRot-y_smooth(1,1)).^2));
            adjacent=sqrt(((T1x_platOriginRot-x_smooth(idist+1,1)).^2)+((T1y_platOriginRot-y_smooth(idist+1,1)).^2));
            trajectoryHeading(i,1)=acosd(adjacent/hypotenuse);
            trajectoryRadians(i,1)=acos(adjacent/hypotenuse);
        end
end

% Compile and export outcome measure
name = tab{j};
if j==1 || j==2 || j==3 || j==4 || j==5;
    %Compile Outcome measures - Trianing
    mwmTrainData=[newC num2cell(D) num2cell(PathLength) num2cell(swimSpeed) num2cell(prefScoreTotal)];
    mwmTrainOutput=cell2table(mwmTrainData,'VariableNames',{'subjectID','Prox_Score','Path_Length','Swim_Speed','Pref_score'});
    %create scatter  
%     scatterGroups=[1 1 1 1 1 1 1 2 2 2 2 2 2 2 1 2 3 3 3 3 3 3 4 4 4 4 4 4 1 1 1 1 1 1 1 2 2 2 2 2 2 2 1 2 3 3 3 3 3 3 4 4 4 4 4 4 1 1 1 1 1 1 1 2 2 2 2 2 2 2 1 2 3 3 3 3 3 3 4 4 4 4 4 4 1 1 1 1 1 1 1 2 2 2 2 2 2 2 1 2 3 3 3 3 3 3 4 4 4 4 4 4 ]';%1=TgM,2=TgF,3=WTM, 4=WTF
%     figure(j); subplot(2,2,1)
%     h=notBoxPlot(D(:),scatterGroups, 'jitter', 0.5,'markMedian',true);
%     title('Search Error','FontSize', 12, 'Color','k','FontName','Helvetica');
% %     set(h(:).mu,'Visible','off')
% %     set(h(:).sd,'Visible','off')
% %     set(h(:).data,'markerfacecolor',)
%     hold on 
%     subplot(2,2,2)
%     notBoxPlot(PathLength(:),scatterGroups, 'jitter', 0.5,'markMedian',true)
%     title('Path Length','FontSize', 12, 'Color','k','FontName','Helvetica');
%     hold on 
%     subplot(2,2,3)
%     notBoxPlot(swimSpeed(:),scatterGroups, 'jitter', 0.5,'markMedian',true)
%     title('Swim Speed','FontSize', 12, 'Color','k','FontName','Helvetica');
%     hold on 
%     subplot(2,2,4)
%     notBoxPlot(prefScoreTotal(:),scatterGroups, 'jitter', 0.5,'markMedian',true)
%     title('Preference Score','FontSize', 12, 'Color','k','FontName','Helvetica');
%     hold off
%     saveas(figure(j),['C:\Users\Ben Clark''s Lab\Google Drive\MWM_Data\Analysis\AnalyzedData\Training\T1\T1',tab{j},'.jpeg']);
    if ismac ==1
        cd '/Users/lauraberkowitz/Google Drive/MWM_Data/Analysis';
    else
        cd 'C:\Users\Ben Clark''s Lab\Google Drive\MWM_Data\Analysis\AnalyzedData\Training';
    end
    writetable(mwmTrainOutput,[pwd filesep 'data_Train' name '.xlsx']);
elseif j==6
    %Compile Outcome measures - Probe Trial
    mwmProbData=[newC num2cell(ProbeD) num2cell(ProbeDf30) num2cell(ProbeDl30) num2cell(ProbepathLength) num2cell(ProbeswimSpeed) num2cell(prefScoreFirst) num2cell(prefScoreLast) num2cell(prefScoreTotal) num2cell(opQuadFirst) num2cell(opQuadLast) num2cell(opQuadTotal) num2cell(timePlatT)...
        num2cell(timePlatA) num2cell(timePlatB) num2cell(timePlatC) num2cell(ProbetrajectoryHeading)];
    mwmProbeOutput=cell2table(mwmProbData,'VariableNames',{'Trial_ID','Prox_Score','Prox_Score_f30','Prox_Score_l30','Path_Length','Swim_Speed','Pref_ScoreFirst', 'Pref_ScoreLast','Pref_ScoreTotal','opQuad_First', 'opQuad_Last','opQuad_Total','Time_PlatLocT'...
        'Time_PlatLocA', 'Time_PlatLocB', 'Time_PlatLocC', 'Trajectory_heading'});
    if ismac ==1
        cd '/Users/lauraberkowitz/Google Drive/MWM_Data/Analysis';
    else
        cd 'C:\Users\Ben Clark''s Lab\Google Drive\MWM_Data\Analysis\AnalyzedData\Probe';
    end
    writetable(mwmProbeOutput,[pwd filesep 'data_Probe' name '.xlsx']);
    saveas(FigTot,'Probe Total Paths','jpeg')
    saveas(Figf30,'Probe first 30s Paths','jpeg')
    saveas( Figl30,'Probe last 30s Paths','jpeg')
    saveas( Figf3,'Probe first 3s Paths','jpeg')
elseif j==7 || j==8 ||j==9 || j==10 || j==11;
    mwmM2PData=[newC num2cell(D) num2cell(Dreference) num2cell(Dprior) num2cell(PathLength) num2cell(swimSpeed) num2cell(trajectoryHeading)];
    mwmM2POutput=cell2table(mwmM2PData,'VariableNames',{'subjectID','Prox_Score','Prox_Ref','Prox_Prior','Path_Length','Swim_Speed','Heading_error'});
    if ismac ==1
        cd '/Users/lauraberkowitz/Google Drive/MWM_Data/Analysis';
    else
        cd 'C:\Users\Ben Clark''s Lab\Google Drive\MWM_Data\Analysis\AnalyzedData\M2P';
    end
    writetable(mwmM2POutput,[pwd filesep 'data_M2P' name '.xlsx']);
elseif j==11
    mwmCuedAData=[newC num2cell(D) num2cell(PathLength) num2cell(swimSpeed) num2cell(trajectoryHeading)];
    mwmCuedAOutput=cell2table(mwmCuedAData,'VariableNames',{'subjectID','Prox_Score','Path_Length','Swim_Speed','Heading_error'});
    if ismac ==1
        cd '/Users/lauraberkowitz/Google Drive/MWM_Data/Analysis';
    else
        cd 'C:\Users\Ben Clark''s Lab\Google Drive\MWM_Data\Analysis\AnalyzedData\T1_Cued';
        writetable(mwmCuedAOutput,[pwd filesep 'data_CuedA' name '.xlsx']);
    end
elseif j==12
    mwmCuedBData=[newC num2cell(D) num2cell(PathLength) num2cell(swimSpeed) num2cell(trajectoryHeading)];
    mwmCuedBOutput=cell2table(mwmCuedBData,'VariableNames',{'subjectID','Prox_Score','Path_Length','Swim_Speed','Heading_error'});
    if ismac ==1
        cd '/Users/lauraberkowitz/Google Drive/MWM_Data/Analysis';
    else
        cd 'C:\Users\Ben Clark''s Lab\Google Drive\MWM_Data\Analysis\AnalyzedData\T1_Cued';
    end
    writetable(mwmCuedBOutput,[pwd filesep 'data_CuedB' name '.xlsx']);
end
clear
end
disp ('Finished with T1 Tracking')

toc


