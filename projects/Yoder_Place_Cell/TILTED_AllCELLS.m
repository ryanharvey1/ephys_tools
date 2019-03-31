% TILTED_AllCELLS 


*NEVER FINISHED DONT USE*


clear;clc;close all
path='/Users/RyanHarvey/Downloads/Place Cells - Tilted Mice';

sessions=dir(path);
sessions={sessions.name};

sessions(contains(sessions,'._'))=[];
sessions(contains(sessions,'.'))=[];

sessions=sessions';

ALLspkfile=[];
c=1;
for i=1:length(sessions)
   

    cd([path,filesep,sessions{i}])
    spkfile=dir('TT*.txt');
    spkfile={spkfile.name}';
    spkfile(contains(spkfile,'._'))=[];
    spkfile(strcmp(spkfile,'.') | strcmp(spkfile,'..'))=[];
    
    for ii=1:length(spkfile)
        ALLspkfile{c}=[path,filesep,sessions{i},filesep,char(spkfile{ii})];
        c=c+1;
    end
    
end
ALLspkfile=ALLspkfile';



control={'SR019','sr019','SR020','sr020'};
tilted={'SR011','SR012','sr012','sr011','RH048','rh048','rh026','RH026','SR021','SR024','SK107'};

controlpath=ALLspkfile(contains(ALLspkfile,control));
tiltedpath=ALLspkfile(contains(ALLspkfile,tilted));

ALLspkfile(~contains(ALLspkfile,control) & ~contains(ALLspkfile,tilted))

[resultsC,ratemapC,dataC,OCCC,sessionC]=meanVELO(Sess1con,ControlSpike,control_paths,'Control',1);
[resultsT,ratemapT,dataT,OCCT,sessionT]=meanVELO(Sess1tilt,TiltedSpike,control_paths,'Tilted',1);

function [results,ratemap,data,OCC,session]=meanVELO(Sess1,spikefile,paths,group,Fig)
FigureLocation=['/Users/RyanHarvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/Figures_For_Paper/',group,'RateMaps'];
ii=1;
for i=1:length(Sess1)
    % Import video file
    split=strsplit(Sess1{i},'ASCII');
    disp(split(1))
    fileID=fopen(char(strcat(split(1),'VT1_DupRem_PrcFields_lights.txt')),'r');
    Intro = textscan(fileID,'%s%s%s');
    fclose(fileID);
    temppath(i)=split(1);
    
    % Import Spikes
    fileID=fopen(spikefile{i},'r');
    if fileID==-1
        disp(['Failed to open',spikefile{i}])
        continue;
    end
    spike = textscan(fileID,'%s');
    fclose(fileID);
    
    % locate TXYS
    T = [Intro{1}];
    X = [Intro{2}];
    Y = [Intro{3}];
    S = [spike{1}];
    
    % Remove ',' & replace empties with NaN
    %     TS=regexprep(T(4:end,1),',','','emptymatch');
    %     VIDEO=regexprep([X(4:end,1),Y(4:end,1)],',','','emptymatch');
    %     empties = cellfun('isempty',VIDEO);
    %     VIDEO(empties) = {NaN};
    
    %     TS=regexprep(TS,'R','','emptymatch');
    %     tempdat=[TS,VIDEO];
    %     tempdat=tempdat(~cellfun('isempty',tempdat(:,1)),:);
    % %     empties = cellfun('isempty',tempdat);
    %     tempdat(cellfun('isempty',tempdat)) = {NaN};
    %     tempdat = sscanf(CStr2String(tempdat, '*'), '%f*');
    %     tempdat = reshape(tempdat,[],3);
    %
    
    tempdat=[regexprep(regexprep(T(4:end,1),',','','emptymatch'),'R','','emptymatch'),regexprep([X(4:end,1),Y(4:end,1)],',','','emptymatch')];
    tempdat=tempdat(~cellfun('isempty',tempdat(:,1)),:);
    tempdat(cellfun('isempty',tempdat))={NaN};
    tempdat=reshape(sscanf(CStr2String(tempdat, '*'), '%f*'),[],3);
    
    %     % remove Rs from TS
    %     TS=regexprep(TS,'R','','emptymatch');
    %     TS=TS(~cellfun('isempty',TS));
    %
    %     % Convert to double
    %     TS = sscanf(CStr2String(TS, '*'), '%f*');
    %     N = sscanf(CStr2String(VIDEO, '*'), '%f*');B = reshape(N,[],2);
    S = sscanf(CStr2String(S, '*'), '%f*');
    
    % REMOVE BETWEEN FRAME NANS WHILE PRESERVING TRACKING ERRORS TO REPRODUCE JENNI'S IGOR CODE
    %     index=[];
    %     for n=2:length(B)
    %         if isnan(B(n))==1 && isnan(B(n-1))==0
    %            index=[index;n];
    %         end
    %     end
    %     B(index,:)=[];
    
    %     interpolate out NaNs from XY
    nanindex=isnan(tempdat(:,2));
    tempdat(isnan(tempdat))=0;
    tempdat(:,5)=nanindex;
    try
        [x,y] = InterpolateTrackerZeros_v2(tempdat(:,2),tempdat(:,3));
    catch
        Z=find(tempdat(:,2)>0);
        Z=(Z(1,1));
        tempdat(1:Z-1,:)=[]; %tempdat(1:Z-1,:)=[];
        [x,y] = InterpolateTrackerZeros_v2(tempdat(:,2),tempdat(:,3));
    end
    
    B=[x,y];
    fail=0;
    while fail==0
        try
            if length(B)>length(tempdat)
                B(end,:)=[];
            elseif length(B)<length(tempdat)
                tempdat(end,:)=[];
            end
            tempdat(:,2:3)=B;
            fail=2;
        catch
            fail=0;
        end
    end
    
    % Restrict path to inside cylinder
    xmax=565; %540;
    xmin=255; %260;
    ymax=430; %405;
    ymin=125; %150;
    x=median([xmax,xmin]); y=median([ymax,ymin]); rad=(min(xmax-xmin,ymax-ymin))/2;
    th = 0:pi/179.5:2*pi; % 0 to 2*pi(6.28318530717959) at 0.0175 increments to equal 360 points
    xunit = rad * cos(th) + x; yunit = rad * sin(th) + y;
    [in,on] = inpolygon(tempdat(:,2),tempdat(:,3),xunit,yunit);
    tempdat(~in,:)=[];
    
    % smooth filtered data
    %     tempdat(:,2:3)=[SmoothVector(tempdat(:,2),30),SmoothVector(tempdat(:,3),30)];
    
    %     % SMOOTHING RAW XY DATA FROM EACH LED
    %     padDsData=[repmat(tempdat(1,:),30,1); tempdat; repmat(tempdat(end,:),30,1)]; %pad ends with last data point
    %     padDsData=[padDsData(:,1),SmoothVector(padDsData(:,2),10),SmoothVector(padDsData(:,3),10)];
    %     tempdat=(padDsData(30+1:end-30,:)); %remove pad
    
    % add vel and filter >70cm/sec
    [vel_cmPerSec,~,~]=InstaVel([tempdat(:,2),tempdat(:,3)],'no',61,60);
    vel_cmPerSec=[vel_cmPerSec;vel_cmPerSec(end,1)];
    tempdat(:,4)=vel_cmPerSec;
    %     tempdat(tempdat(:,4)>100,:)=[];
    
    TS_new =interp1(tempdat(:,1),tempdat(:,1),S,'linear');
    X=interp1(tempdat(:,1), tempdat(:,2), S, 'linear');
    Y=interp1(tempdat(:,1), tempdat(:,3), S, 'linear');
    V=interp1(tempdat(:,1), tempdat(:,4), S, 'linear');
    NANI=interp1(tempdat(:,1), tempdat(:,5), S, 'linear');
    
    
    % CONCAT, SORT, & ADD SPIKE BINARY
    TSXY=sortrows([[TS_new X Y V NANI ones(size(TS_new,1),1)];[tempdat,zeros(length(tempdat),1)]],1);
    
    % remove spikes with missing xy coordinates
    TSXY(TSXY(:,5)>0,:)=[];
    TSXY(:,5)=TSXY(:,6);
    TSXY(:,6)=[];
    
    % calculate instantaneous firing rate
    spikeVal=TSXY(TSXY(:,5)==1,1);
    
    ifr=IFR(spikeVal,TSXY(:,5),round(length(tempdat)/60),60);
    TSXY(:,6)=ifr;
    
    [ angle ] = XYangle(TSXY(:,2),TSXY(:,3));
    angle=[angle;angle(end,1)];
    TSXY=[TSXY,angle];
    
    TSXY(isnan(TSXY(:,1)),:)=[];
    
    % find start and end times
    TSdiff=diff(TSXY(:,1));
    mul=1;
    while true
        SE=find(TSdiff>std(TSdiff)*mul);
        mul=mul+1;
        if length(SE)==4; break;end
    end
    Start=[1;SE(1)+1;SE(2)+1;SE(3)+1;SE(4)+1];
    End=[SE(1);SE(2);SE(3);SE(4);length(TSXY)];
    
    
    % //////////////// LOOP THROUGH EACH SESSION ////////////////
    sub=[1,3,5,7,9];
    sub2=[2,4,6,8,10];
    subs=1;
    for j=1:5
        Session=TSXY(Start(j):End(j),:);
        % speed and ifr corr
        try [c,p]=corr(Session(Session(:,5)==1,4),Session(Session(:,5)==1,6));catch;c=NaN;p=NaN;end
        
        % collect occ
        if i==1;pass=1;else;if strcmp(temppath(i-1),split(1))==1;pass=0;else;pass=1;end;end
        if j==1 && pass==1;OCC.(['Cell',num2str(ii),'Session',num2str(j)])=Session(:,2:3);end
        
        % calc stats
        [SmoothRateMap,nBinsx,nBinsy,occ,Coherence]=bindataYoder(Session,60,Session(Session(:,5)==1,:),'no',61);
        rY = reshape(SmoothRateMap,nBinsx*nBinsy,1); % reshape data into column
        rY(isnan(rY))=0; rY(isinf(rY))=0;
        occRSHP = reshape(occ,nBinsx*nBinsy,1); % reshape data into column
        occSUM = sum(occRSHP); % summed occupancy
        pX = occRSHP./occSUM; % normalized occupancy
        InformationContent=InformationPerSpike(rY,pX); % from NSMA toolbox
        PeakRate = max(rY);
        OverallFR = mean(rY);
        nSpikes=sum(Session(:,5));
        sparsity = Sparsity(rY'); % from NSMA toolbox
        NumbActiveBins = numel(find(rY > PeakRate*0.20));
        
        SmoothRateMap2=SmoothRateMap;
        [field,FieldWidth]=FindFF2D(SmoothRateMap2);
        FieldWidth=FieldWidth*(61/length(field));
        if isnan(field)==0
            nanPlacement=isnan(SmoothRateMap2);
            SmoothRateMap2(~field)=0; SmoothRateMap2(nanPlacement)=NaN;
            
            Map4bordertemp=SmoothRateMap2;Map4bordertemp(~field)=NaN;
            infieldFR=nanmean(reshape(Map4bordertemp,[],1));
            
            Map4bordertemp=SmoothRateMap;Map4bordertemp(logical(field))=NaN;
            outfieldFR=nanmean(reshape(Map4bordertemp,[],1));
        else
            infieldFR=0;
            outfieldFR=0;
        end
        % find field
        [row,col]=find(field);
        %find boundary
        k=boundary(row,col);
        bound(i).control=[col(k),row(k)];
        E=Eccentricity(col(k),row(k));
        try
            [borderScore,Field2Wall]=BorderScore(SmoothRateMap2);
            Field2Wall=Field2Wall*(61/length(SmoothRateMap2));
        catch
            borderScore=NaN;
            Field2Wall=NaN;
        end
        
        results(ii,:,j)=[PeakRate,nSpikes,OverallFR,NumbActiveBins,sparsity,InformationContent,Coherence,Field2Wall,borderScore,FieldWidth,infieldFR,outfieldFR,E,c,p];
        ratemap.(['Cell',num2str(ii),'Session',num2str(j)])=SmoothRateMap;
        data.(['Cell',num2str(ii),'Session',num2str(j)])=Session;
        
        if Fig==1
            % plot
            fig=figure(ii); fig.Color=[1 1 1];
            subplot(5,2,sub(subs));
            h=plot(Session(:,2), Session(:,3), 'LineWidth', 1, 'color', 'k');
            hold on; axis off
            spks_VEL=Session(Session(:,5)==1,:);
            scatter(spks_VEL(:,2), spks_VEL(:,3), 5, 'filled', 'r');
            box off; axis image
            set(gca, {'YDir'}, {'reverse'});
            SmoothRateMap(isnan(SmoothRateMap))=0;
            upsamRateMap=PerfectCircRateMap(SmoothRateMap,0);
            subplot(5,2,sub2(subs));
            imAlpha=ones(size(upsamRateMap));
            imAlpha(isnan(upsamRateMap))=0;
            imagesc(upsamRateMap,'AlphaData',imAlpha);
            set(gca,'color',[1 1 1]);
            box off
            axis image
            axis off
            colormap jet
            shading interp
            subs=subs+1;
        end
    end
    tempsess=strsplit(split{1},'/');
    session{ii,1}=tempsess{6};
    if Fig==1
        if sum(get(0,'screensize')==[1,1,1440,900])==4
            set(figure(ii),'Position',[1 5 720 800]);
        else
            set(figure(ii),'Position',[1 1 960 990]);
        end
        pause(1)
        print(figure(ii),'-dpng', '-r600',[FigureLocation,filesep,'Cell',num2str(ii),'.png'])
        close all
    end
    ii=ii+1;
end
end





% clear;clc;close all
% path='/Volumes/Ryan_4TB/Place_Cell_Data/Place Cells - Tilted Mice';
% 
% sessions=dir(path);
% sessions={sessions.name};
% 
% sessions(contains(sessions,'._'))=[];
% sessions(contains(sessions,'.'))=[];
% 
% sessions=sessions';
% 
% control={'SR019','sr019','SR020','sr020'};
% tilted={'SR011','SR012','sr012','sr011','RH048','rh048','rh026','RH026','SR021','SR024','SK107'};
% 
% for i=1:length(sessions)
%     files=dir([path,filesep,sessions{i}])
%     files={files.name};
% 
%     files(contains(files,'._'))=[];
%     
%     files(strcmp(files,'.') | strcmp(files,'..'))=[];
%     
%     
% end




