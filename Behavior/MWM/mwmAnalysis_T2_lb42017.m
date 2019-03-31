% WaterMazeAnalysis analyzes watermaze data
% Input: Compiled XY coordinates from FIJI manual track plugin
%      Coordinates should be stored in an excel file named
%      T#_tracking i.e. T1_tracking (Timepoint 1 tracking file)
%
% Created by BC 2014, adapted for FIJI tracking by LB 2016-2017
clc;
close all;
clear;

addpath(genpath('/Users/lauraberkowitz/Google Drive/MATLAB/chronux'))
addpath(genpath('/Users/ryanharvey/GoogleDrive/MatlabDir/chronux_2_11'))
addpath(genpath('C:\Users\Ben Clark''s Lab\Google Drive\MATLAB\chronux'))
% addpath(genpath('/Users/ryanharvey/GoogleDrive/MatlabDir/BClarkToolbox/Analyses/spikeCode/MEX'))
addpath(genpath('C:\Users\Ben Clark''s Lab\Documents\MATLAB\MEX'))

tic
%Import Data
% if ismac == 1
%     trackData=readtable('/Users/lauraberkowitz/Google Drive/MWM_Data/Analysis/T2_Tracking.xlsx');
% else
%     trackData=readtable('F:\Tracking\Analysis\T2_Tracking.xlsx'); %input path to data
% end
trackData=readtable('/Users/ryanharvey/GoogleDrive/MatlabDir/BClarkToolbox/Analysis/MWM/T2_Tracking.xlsx');
cd '/Users/ryanharvey/GoogleDrive/MatlabDir/BClarkToolbox/Analysis/MWM'
%%____________________________Initialize Pool Parameters_____________________

%Pool Parameters
%Pool Origin - obtained using FIJI Measure function
%____T2__________
T2x_Origin=292.5;
T2y_Origin=238.52;

%Platform Origin - polygon function to outline platform, measure to obtain
%origin. Deduct
%____T2__________
T2x_platOrigin=371.497;
T2y_platOrigin=328.283;

%Pool boundaries - x/y min and max must intersect origin - obtained using
%elipse tool (measure tool to get centroid or origin).
%____T2__________
T2x_max=552;
T2x_min=35;
T2y_max=479;
T2y_min=0;

%%##### Convert Coordinates from pixels to cm ############################

%Pixel to cm converstion factor - 150 indicates pool diameter -
%____T2__________
T2x_conv=150/(T2x_max-T2x_min); %cm/pixel
T2y_conv=150/(T2y_max-T2y_min);

%Pool Boundaries in cm
%____T2__________
T2x_max_cm=(T2x_max*T2x_conv);
T2x_min_cm=(T2x_min*T2x_conv);
T2y_max_cm=(T2y_max*T2y_conv);
T2y_min_cm=(T2y_min*T2y_conv);

%Pool Origin in cm
T2x_Origin_cm=T2x_Origin*T2x_conv;
T2y_Origin_cm=T2y_Origin*T2y_conv;

%Platform Origin in cm
T2x_platOrigin_cm=T2x_platOrigin*T2x_conv;
T2y_platOrigin_cm=T2y_platOrigin*T2y_conv;

%Platform Origin transposed relative to pool origin of (0,0)
T2x_platOriginTranspose=T2x_platOrigin_cm-T2x_Origin_cm;
T2y_platOriginTranspose=T2y_platOrigin_cm-T2y_Origin_cm;

%Create squre platform (make lower left corner position
xpos=T2x_platOriginTranspose-8; %for a 16x16cm square platform
ypos=T2y_platOriginTranspose-8;

%Convert tracking data to cm
trackData{:,3}= (((trackData{:,3}))-(T2x_Origin))*T2x_conv;
trackData{:,4}=((trackData{:,4}-(T2y_Origin))*T2y_conv);

%=-=-=-=-=-===-=-=-=--== Rotate data points=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
plat=[T2x_platOriginTranspose;T2y_platOriginTranspose];
% choose a point which will be the center of rotation
x_center = 0;
y_center = 0;
% create a matrix which will be used later in calculations
centerPlat =[x_center; y_center];
centerData = repmat([x_center; y_center], 1, length(trackData{:,3}));
% define a counter-clockwise rotation matrix
theta = deg2rad(-5.8070);     % dervied offline x=atan(o/a)
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
T2x_platOriginRot= rotPlat(1,:);
T2y_platOriginRot= rotPlat(2,:);
%-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-==-Create Platform=-=-=-=-=-=-=-=-=-=-=-=-
%Create squre platform pos for rotated platform (make lower left corner position
xposRot=T2x_platOriginRot-8; %for a 16x16cm square platform
yposRot=T2y_platOriginRot-8;

%% Initialize variables to loop through xy coordinate data
% nameSplit=strsplit((name),'_'); %Identify timepoint via filename
[C, ia]=unique(trackData(:,1),'stable');
newC=reshape(C{:,1},[],1);
% filepath ='F:\Tracking\Analysis\AnalyzedData';

%Figure Parameters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
titles={'Subject 1,male,Tg+','Subject 2,male,Tg+','Subject 3,male,Tg+','Subject 4,male,Tg+',...
    'Subject 5,male,Tg+',...
    'Subject 6,male,Tg+','Subject 7,male,Tg+','Subject 8,female,Tg+','Subject 9,female,Tg+','Subject 10,female,Tg+',...
    'Subject 11,female,Tg+','Subject 12,female,Tg+','Subject 13,female,Tg+','Subject 14,female,Tg+','Subject 15,male,Tg+',...
    'Subject 16,female,Tg+','Subject 17,male,WT','Subject 18,male,WT','Subject 19,male,WT','Subject 20,male,WT',...
    'Subject 21,male,WT','Subject 22,male,WT','Subject 23,female,WT','Subject 24,female,WT','Subject 25,female,WT',...
    'Subject 26,female,WT','Subject 27,female,WT','Subject 28,female,WT'};
%Pool Parameters
x=0;
y=0;
r=165/2; %150 is pool diameter
ang=0:0.01:2*pi;
xp=r*cos(ang);
yp=r*sin(ang);

%Create figure
FigTot=figure('Name','Time Point 2 Probe Trial','NumberTitle','off'); %change title to fit data description
Figf30=figure('Name','Time Point 2 Probe Trial: First 30s','NumberTitle','off'); %change title to fit data description
Figl30=figure('Name','Time Point 2 Probe Trial: Last 30s','NumberTitle','off');%change title to fit data description
Figf3=figure('Name','Time Point 2 Probe Trial: Initial Trajectory','NumberTitle','off');
%Plot Smoothed Data~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Initilize subgroup matrices for heat maps/scatter box plots
smoothTgM=[];
smoothTgF=[];
smoothWTM=[];
smoothWTF=[];
smoothTg=[];
smoothWT=[];

nBins = 65; % adjusted so the bins are ~5x5cm as in Koenig et al 2011
MinY = min(trackData{:,4});
MaxY = max(trackData{:,4});
MinX = min(trackData{:,3});
MaxX = max(trackData{:,3});
edges{1} = linspace(MinY, MaxY, nBins+1);
edges{2} = linspace(MinX, MaxX, nBins+1);

%Loop through each trial, process/plot data and compute output variables.
for i=1:size(C,1) % %loop through #of participants
    if i==29; break; end
    try
        trialDsnuts=trackData(ia(i):(ia((i+1))-1),:);%Use unique index to process one session at a time.
    catch
    end
    if i==length(ia);
        trialDsnuts=trackData(ia(i):end,:);
    end
    trialData=reshape(trialDsnuts{:,2:4},[],3);%convert to double for processing with resample/runline
    %     trialDataPad=[zeros(5,3); trialData; zeros(5,3)];
    %     Change sampling rate to 10Hz
    %     -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-CHANGE PER TIMEPOINT
    T2Fs=11; %observed frame rate in Hz
    T2dsFs=10; %desired frame rate in Hz
    [p,q]=rat(T2dsFs/T2Fs);
    dsData=resample(trialData,p,q);
    
    %     dsDataDown=downsample(trialData,2);
    trial=(strsplit(newC{i,:},'_'));
    if strcmp('t21',trial{2})%Probe session - shrink to 60sec duration
        dsData((length(dsData)-(((length(dsData)/10)-60)*10))+1:end,:)=[];
    end
    %Smooth Data
    x_smooth=runline(dsData(2:end,2),3,1); %Smooth data
    y_smooth=runline(dsData(2:end,3),3,1); %Smooth data
    
    % STRATEGY ANALYSIS RYAN HARVEY 
    % exports csv for every session

    Samples=1:size(x_smooth,1); TS=0:1/10:size(x_smooth,1);    TS=TS(1:size(x_smooth,1));
    id=table2array(trialDsnuts(1,1)); id=strsplit(id{1},'_');
    subject=strsplit(id{3},'j'); subject=subject(2);
    trialInfo={'Time_Point',num2str(id{1}),NaN,NaN;'Sample',num2str(size(x_smooth,1)),NaN,NaN;'ID',num2str(subject{1}),NaN,NaN;'trial',id{2},NaN,NaN};
        
    array_export=num2cell([Samples',TS',x_smooth,y_smooth]);
    
    exportsheet=[trialInfo;{'Samples','TS','X','Y'};array_export];
    
    filename=[pwd filesep 'Time_Point' num2str(id{1}) filesep strjoin({'day1','_00',num2str(i)},'') '.csv'];
    
    if i==1; mkdir([pwd filesep 'Time_Point' num2str(id{1})]); end

    fid=fopen(filename,'w');
    fprintf(fid,'%s, %s, %s, %s\n',exportsheet{1,:});
    fprintf(fid,'%s, %s, %s, %s\n',exportsheet{2,:});
    fprintf(fid,'%s, %s, %s, %s\n',exportsheet{3,:});
    fprintf(fid,'%s, %s, %s, %s\n',exportsheet{4,:});
    fprintf(fid,'%s, %s, %s, %s\n',exportsheet{5,:});
    fclose(fid);
    dlmwrite(filename,exportsheet(6:end,:),'-append');
    
    
    %-=-=-=-=-=-=-=-=-=-=-=-=-=-=-Obtain Outcome Measures Training-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    %Directional heading, path length, swim speed, circular statistics
    %for trials 9 - 20 (training)
    %Directional heading
    
    if strcmp('t21',trial{2})==0;
        %Calculate initial trajectory
%         Locate xy up to 45cm
%         repTraject=repmat([x_smooth(1,1) y_smooth(1,1)],([(size(x_smooth,1)-1) 1]));
%         trajectoryIndex=[x_smooth(2:end,1) y_smooth(2:end,1)]- repTraject(:,:); %compare each point to start
%         trajectoryDistance=sqrt((trajectoryIndex(:,1).^2)+(trajectoryIndex(:,2).^2)); %calculate distance
%         for idist=1:length(trajectoryDistance)
%             if trajectoryDistance(idist)>=44.5; break; end; %break when the difference is at least 44.5cm
%         end
%         hypotenuse=sqrt(((T2x_platOriginTranspose-x_smooth(1,1)).^2)+((T2y_platOriginTranspose-y_smooth(1,1)).^2));
%         adjacent=sqrt(((T2x_platOriginTranspose-x_smooth(idist,1)).^2)+((T2y_platOriginTranspose-y_smooth(idist,1)).^2));
%         trainingTrajectoryRadians(i,1)=real(acos(adjacent/hypotenuse)); %in radians for circular statistics
%         trainingTrajectoryHeading(i,1)=real(acosd(adjacent/hypotenuse));
%         cosRatio=(adjacent/hypotenuse);

        %Calculate average proximity
        repTest=repmat([T2x_platOriginRot T2y_platOriginRot],([size(x_smooth,1) 1])); %initialize matrix for platform location
        proximity=[x_smooth y_smooth]- repTest(:,:); %calculate difference between path and platform location
        TrainingD(i,1)= (sum(sqrt((proximity(:,1).^2)+(proximity(:,2).^2))))-(sqrt(((x_smooth(1,1)-(T2x_platOriginRot))^2)+...
            (y_smooth(1,1)-(T2y_platOriginRot))^2)); %cumulative Proximity sans distance from start to plat origin. 
        
        
        %Calculate Path length & Swim Speed
        sqrXDiff=(diff(x_smooth)).^2;
        sqrYDiff=(diff(y_smooth)).^2;
        pathDist=sqrt(sqrXDiff + sqrYDiff);
        trainingPathLength(i,1)=sum(pathDist); %Total Path length
        trainingswimSpeed(i,1)=(trainingPathLength(i,1)/(size(x_smooth,1)/10)); %total swim speed
        
        %Calculate preference score total
        timeC=(sum(x_smooth(:)<=0 & y_smooth(:)>=0))/10; %calculate data points in target quadrant in seconds
        timeT=(sum(x_smooth(:)>=0 & y_smooth(:)>=0))/10;
        timeA=(sum(x_smooth(:)>=0 & y_smooth(:)<=0))/10;
        timeB=(sum(x_smooth(:)<=0 & y_smooth(:)<=0))/10;
        prefScoreTrainTotal(i,1)=((timeT-timeA)+(timeT-timeB)+(timeT-timeC)/3);
        
    elseif strcmp('t21',trial{2})==1;
        
        %Plot data, add to figure whole session
        set(0,'CurrentFigure',FigTot);subplot(4,7,i)
        plot(x+xp,y+yp);
        %     axis([-150 150]);
        title(titles(i),'FontSize', 12, 'Color','k','FontName','Helvetica');
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
        title(titles(i),'FontSize', 12, 'Color','k','FontName','Helvetica');
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
        title(titles(i),'FontSize', 12, 'Color','k','FontName','Helvetica');
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
        title(titles(i),'FontSize', 12, 'Color','k','FontName','Helvetica');
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
        sbj=regexp(trial{3},'\d*', 'match');
        TgM={'1' '2' '3' '4' '5' '6' '7' '15'}; %subject IDs
        TgF={'8' '9' '10' '11' '12' '13' '14' '16'};
        WTM={'17' '18' '19' '20' '21' '22'};
        WTF={'23' '24' '25' '26' '27' '28'};
        Tg={'1' '2' '3' '4' '5' '6' '7' '15' '8' '9' '10' '11' '12' '13' '14' '16'};
        WT={'17' '18' '19' '20' '21' '22' '23' '24' '25' '26' '27' '28'};
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
        
        %Compile data for heat maps between groups - first 30s i.e. first 300
        %rows.
        if ismember(sbj,Tg);
            %         catSubTgM = repmat(sbj,length(x_smooth),1);
            iTg= [x_smooth(1:300,:) y_smooth(1:300,:)]; %concatenate first 30s
            smoothTg=[smoothTg; iTg];
        else ismember(sbj,WT);
            %         catSubTgF = repmat(sbj,length(x_smooth),1);
            iWT=[x_smooth(1:300,:) y_smooth(1:300,:)];%concatenate first 30s
            smoothWT=[smoothWT; iWT];
        end
        
        %~~~~~~~~~~~~~~~~~~~  Calculate Outcome Measures  ~~~~~~~~~~~~~~~~~~~~~~~
        %Calculate average proximity
        repTest=repmat([T2x_platOriginRot T2y_platOriginRot],([size(x_smooth,1) 1])); %initialize matrix for platform location
        proximity=[x_smooth y_smooth]- repTest(:,:); %calculate difference between path and platform location
        D(i,1)= mean(sqrt((proximity(:,1).^2)+(proximity(:,2).^2))); %Average Proximity entire trial
        Df30(i,1)=mean(sqrt((proximity(1:300,1).^2)+(proximity(1:300,2).^2)));
        Dl30(i,1)=mean(sqrt((proximity(301:end,1).^2)+(proximity(301:end,2).^2)));
        
        %Create path differences between points
        sqrXDiff=(diff(x_smooth)).^2;
        sqrYDiff=(diff(y_smooth)).^2;
        pathDist=sqrt(sqrXDiff + sqrYDiff);
        pathLength(i,1)=sum(pathDist); %Total Path length
        swimSpeed(i,1)=(pathLength(i,1)/(size(x_smooth,1)/10));
        
        %Calculate initial trajectory
        % Locate xy up to 45cm
        repTraject=repmat([x_smooth(1,1) y_smooth(1,1)],([(size(x_smooth,1)-1) 1]));
        trajectoryIndex=[x_smooth(2:end,1) y_smooth(2:end,1)]- repTraject(:,:);
        trajectoryDistance=sqrt((trajectoryIndex(:,1).^2)+(trajectoryIndex(:,2).^2));
        for idist=1:length(trajectoryDistance)
            if trajectoryDistance(idist)>=44.5; break; end;
            %        distance(idist)=trajectoryIndex(idist);
            %        if distance>=44.5; break; end
        end
        hypotenuse=sqrt(((T2x_platOriginRot-x_smooth(1,1)).^2)+((T2y_platOriginRot-y_smooth(1,1)).^2));
        adjacent=sqrt(((T2x_platOriginRot-x_smooth(idist+1,1)).^2)+((T2y_platOriginRot-y_smooth(idist+1,1)).^2));
        trajectoryHeading(i,1)=acosd(adjacent/hypotenuse);
        trajectoryRadians(i,1)=acos(adjacent/hypotenuse);
        
        %     % Index xy to 45cm
        %     x4angle=x_smooth(idist+1,1); y4angle=y_smooth(1:idist+1,1);
        %     % Extract angle up to 20cm
        %     for iangle=1:length(x4angle)
        %         if (iangle+1)>length(x4angle); break; end
        %         ExtractedAngle(iangle) = atan2d(y4angle(iangle+1)-y4angle(1),x4angle(iangle+1)-x4angle(1)) + 360*((y4angle(iangle+1)-y4angle(1)<0));
        %         %     ExtractedAngle = atan2(table(:,5)-table(:,3),table(:,4)-table(:,2)) + 360*((table(:,5)-table(:,3))<0);
        %     end
        %     directTraject(i,1)=atan2d(y4angle(1)-yposRot,x4angle(1)-xposRot) + 360*((y4angle(1)-yposRot)<0);
        %     meanTrajectory(i,1)=mean(ExtractedAngle);
        %         clear repTraject trajectoryIndex trajectoryDistance idist Hypotenuse adjacent
        
        %Calculate time in platform area - Quadrant T
        timePlatT(i,1)=(sum(x_smooth>=(T2x_platOriginRot-8) & x_smooth<=(T2x_platOriginRot+8) &...
            y_smooth>=(T2y_platOriginRot-8) & y_smooth<=(T2y_platOriginRot+8)))/10;
        
        %Calculate time in platform area - Quadrant A
        timePlatA(i,1)=(sum(x_smooth>=((T2x_platOriginRot))-8 & x_smooth<=((-1.*T2x_platOriginRot)+8) &...
            y_smooth>=((-1.*T2y_platOriginRot)-8) & y_smooth<=((-1.*T2y_platOriginRot)+8)))/10;
        
        %Calculate time in platform area - Quadrant B
        timePlatB(i,1)=(sum(x_smooth>=((-1.*(T2x_platOriginRot))-8) & x_smooth<=((-1.*(T2x_platOriginRot))+8) &...
            y_smooth>=((-1.*(T2y_platOriginRot))-8) & y_smooth<=((-1.*T2y_platOriginRot)+8)))/10;
        
        %Calculate time in platform area - Quadrant C
        timePlatC(i,1)=(sum(x_smooth>=((-1.*(T2x_platOriginRot))-8) & x_smooth<=((-1.*(T2x_platOriginRot))+8) &...
            y_smooth>=(T2y_platOriginRot-8) & y_smooth<=(T2y_platOriginRot+8)))/10;
        
        
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
    end
end
%Heat maps
nBins = 64; % adjusted so the bins are ~5x5cm as in Koenig et al 2011
MinY = min(trackData{:,4});
MaxY = max(trackData{:,4});
MinX = min(trackData{:,3});
MaxX = max(trackData{:,3});
edges{1} = linspace(MinY, MaxY, nBins+1);
edges{2} = linspace(MinX, MaxX, nBins+1);

figure(5)
subplot(1,4,1);
binTgM=hist3(smoothTgM,'Edges',edges);
binTgM(1,:) = [];
binTgM(:,end) = [];
[ circBinMapTgM ] = circRateMap(binTgM,nBins );
h=pcolor(circBinMapTgM);
hold on
colormap(hot)%flipud to reverse colormap gradient
axis off
title('Tg+ Males' ,'FontSize', 20, 'Color','k','FontName','Times New Roman');
box off
set(h,'EdgeColor','none');
axis image
% set(gca,'YDir','reverse')
hold off

figure (5)
b=subplot(1,4,2);
binTgF=hist3(smoothTgF,'Edges',edges);
binTgF(1,:) = [];
binTgF(:,end) = [];
[ circBinMapTgF ] = circRateMap(binTgF,nBins );
h=pcolor(circBinMapTgF);
colormap(hot)%flipud to reverse colormap gradient
title('Tg+ Females' ,'FontSize', 20, 'Color','k','FontName','Times New Roman');
shading interp
axis off
box off
set(h,'EdgeColor','none');
axis image
% set(gca,'YDir','reverse')
hold off

figure (5)
c=subplot(1,4,3);
binWT=hist3(smoothWTM,'Edges',edges);
binWT(1,:) = [];
binWT(:,end) = [];
[ circBinMapWTM ] = circRateMap(binWT,nBins );
h=pcolor(circBinMapWTM);
colormap(hot)%flipud to reverse colormap gradient
title('WT Males' ,'FontSize', 20, 'Color','k','FontName','Times New Roman');
shading interp
axis off
box off
set(h,'EdgeColor','none');
axis image
% set(gca,'YDir','reverse')
hold off

figure (5)
d=subplot(1,4,4);
binWTF=hist3(smoothWTF,'Edges',edges);
binWTF(1,:) = [];
binWTF(:,end) = [];
[ circBinMapWTF ] = circRateMap(binWTF,nBins );
h=pcolor(circBinMapWTF);
colormap(hot)%flipud to reverse colormap gradient
title('WT Females' ,'FontSize', 20, 'Color','k','FontName','Times New Roman');
shading interp
axis off
box off
set(h,'EdgeColor','none');
axis image
% set(gca,'YDir','reverse')
hold off

figure(6)
%Plot heatmaps
d=subplot(1,2,1);
binTg=hist3(smoothTg,'Edges',edges);
binTg(1,:) = [];
binTg(:,end) = [];
[ circBinMapTg ] = circRateMap(binTg,nBins );
h=pcolor(circBinMapTg);
colormap(hot)%flipud to reverse colormap gradient
title('Tg+' ,'FontSize', 20, 'Color','k','FontName','Times New Roman');
shading interp
axis off
box off
set(h,'EdgeColor','none');
axis image
% set(gca,'YDir','reverse')
hold off

%Plot heatmaps
figure(6)
d=subplot(1,2,2);
binWT=hist3(smoothWT,'Edges',edges);
binWT(1,:) = [];
binWT(:,end) = [];
[ circBinMapWT ] = circRateMap(binWT,nBins );
h=pcolor(circBinMapWT);
colormap(hot)%flipud to reverse colormap gradient
title('WT' ,'FontSize', 20, 'Color','k','FontName','Times New Roman');
shading interp
axis off
box off
set(h,'EdgeColor','none');
axis image
% set(gca,'YDir','reverse')
hold off

%compute mean head angle for groups
[TgMale_meanHeading,TgMale_ul,TgMale_ll]=circ_mean(trajectoryRadians(1:7 & 15,1),[],1);
[TgFemale_meanHeading,TgFemale_ul,TgFemale_ll]=circ_mean(trajectoryRadians(8:14 & 16,1),[],1);
[WTM_meanHeading,WTM_ul,WTM_ll]=circ_mean(trajectoryRadians(17:22,1),[],1);
[WTF_meanHeading,WTF_ul,WTF_ll]=circ_mean(trajectoryRadians(23:28,1),[],1);

%% Compile and export outcome measures
%Compile Outcome measures - Trianing
mwmTrainData2Cell=[table2cell(C) num2cell(TrainingD) num2cell(trainingPathLength) num2cell(trainingswimSpeed) num2cell(prefScoreTrainTotal)];
mwmTrainData=[mwmTrainData2Cell(29:end,1) mwmTrainData2Cell(29:end,2) mwmTrainData2Cell(29:end,3) mwmTrainData2Cell(29:end,4) mwmTrainData2Cell(29:end,5)];
mwmTrainOutput=cell2table(mwmTrainData,'VariableNames',{'subjectID','Prox_Score','Path_Length','Swim_Speed','Pref_score'});

%Compile Outcome measures - Probe Trial
[path,name]=fileparts('F:\Tracking\Analysis\T2_Tracking.xlsx'); %Get name of current file
baseFileName = name;
% mwmSmoothData = [smoothTgM smoothTgF smoothWTM smoothWTF];
% smoothOutput=cell2table(mwmSmoothData,'Variable Names',{'Tg+_M_SubID','x_smooth','y_smooth','Tg+_F_SubID','x_smooth','y_smooth','WT_M_SubID','x_smooth','y_smooth','WT_F_SubID','x_smooth','y_smooth'});
mwmProbData=[table2cell(C(1:28,1)) num2cell(D) num2cell(Df30) num2cell(Dl30) num2cell(pathLength) num2cell(swimSpeed) num2cell(prefScoreFirst) num2cell(prefScoreLast) num2cell(prefScoreTotal) num2cell(opQuadFirst) num2cell(opQuadLast) num2cell(opQuadTotal) num2cell(timePlatT)...
    num2cell(timePlatA) num2cell(timePlatB) num2cell(timePlatC) num2cell(trajectoryHeading)];
mwmProbeOutput=cell2table(mwmProbData,'VariableNames',{'Trial_ID','Prox_Score','Prox_Score_f30','Prox_Score_l30','Path_Length','Swim_Speed','Pref_ScoreFirst', 'Pref_ScoreLast','Pref_ScoreTotal','opQuad_First', 'opQuad_Last','opQuad_Total','Time_PlatLocT'...
    'Time_PlatLocA', 'Time_PlatLocB', 'Time_PlatLocC', 'Trajectory_heading'});
if ismac ==1
    cd '/Users/lauraberkowitz/Google Drive/MWM_Data/Analysis';
else
    cd 'F:\Tracking\Analysis\AnalyzedData';
end

%Compile angle data

writetable(mwmProbeOutput,[pwd filesep 'data_Probe' name '.xlsx']);
writetable(mwmTrainOutput,[pwd filesep 'data_Train' name '.xlsx']);
% writetable(smoothOutput,['smoothedCoords' name '.xlsx']);
fprintf(1, 'Finished with %s\n', baseFileName);
toc






