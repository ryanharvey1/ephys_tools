% LOAD_NTT_FILES
% MeanRateMaps
clear ; close all ;
addpath(genpath('/Users/RyanHarvey/GoogleDrive/MatlabDir/NSMA_Toolbox'));
addpath(genpath('/Users/RyanHarvey/GoogleDrive/MatlabDir/BClarkToolbox/Analysis/Yoder_Place_Cell'));
addpath('/Users/RyanHarvey/GoogleDrive/MatlabDir/BClarkToolbox/Analysis');
addpath('/Users/RyanHarvey/GoogleDrive/MatlabDir/CircStat2012a')
addpath('/Users/ryanharvey/GoogleDrive/MatlabDir/CStr2String')
addpath('/Users/RyanHarvey/GoogleDrive/MatlabDir/BayesRR')
% FigureLocation='/Users/RyanHarvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/Figures_For_Paper';
% paths=importdata('/Users/RyanHarvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/GW- Field_StatsPlace Cells_Tilted_Mice_1_1.xlsx');
load('/Users/RyanHarvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/GW_Field_StatsPlaceCells_Tilted_Mice.mat')
% load('/Users/RyanHarvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/NewData.mat')
path_group=data.textdata.Field_StatsPlaceCells0x2DTilted(:,1:2);
path_group(1,:)=[];
group=path_group(:,2);
control_paths=path_group(strcmp(group,'Control'),1);
tilted_paths=path_group(strcmp(group,'Tilted'),1);


% just first session
Sess1con=control_paths(1:5:length(control_paths));
Sess1tilt=tilted_paths(1:5:length(tilted_paths));

% FOR SPIKE FILES
ControlSpike=regexprep(Sess1con,'.txtrmap1.txt','.txt');
TiltedSpike=regexprep(Sess1tilt,'.txtrmap1.txt','.txt');

% % FOR RATEMAPS
% Controlratemap=regexprep(control_paths,'_Data.mat','.txt');
% Tiltedratemap=regexprep(tilted_paths,'_Data.mat','.txt');

% [PeakRate,nSpikes,OverallFR,NumbActiveBins,sparsity,InformationContent,Coherence,Field2Wall,borderScore,FieldWidth,infieldFR,outfieldFR,E,c,p];


[resultsC,ratemapC,dataC,OCCC,sessionC]=meanVELO(Sess1con,ControlSpike,control_paths,'Control',0);
[resultsT,ratemapT,dataT,OCCT,sessionT]=meanVELO(Sess1tilt,TiltedSpike,control_paths,'Tilted',0);

results.control.results=resultsC;
results.control.ratemap=ratemapC;
results.control.data=dataC;
results.control.occ=OCCC;
results.control.id=sessionC;

results.tilted.results=resultsT;
results.tilted.ratemap=ratemapT;
results.tilted.data=dataT;
results.tilted.occ=OCCT;
results.tilted.id=sessionT;

% % direction
% done=PlotsStat(HDc.mean_vector_length(:,:,:),HDt.mean_vector_length(:,:,:),{'MeanVectorLength'})
%
% [ AllStats ] = ScatterBox([HDc.mean_vector_length(:,:,1),HDc.Direct_infoContent(:,:,1),HDc.preferred_Direction(:,:,1)],[HDt.mean_vector_length(:,:,1),HDt.Direct_infoContent(:,:,1),HDt.preferred_Direction(:,:,1)],{'Control','Tilted'},{'MeanVectorLength','DirectInfoContent','PreferredDirection'},2)

% Polar_Fig=figure;
% Polar_Fig.Color=[1 1 1];
% p1=polarhistogram(deg2rad(HDc.preferred_Direction(:,:,1)),60,'FaceColor','k','FaceAlpha',0.8,'EdgeColor','w','Normalization','probability'); hold on
% p2=polarhistogram(deg2rad(HDt.preferred_Direction(:,:,1)),60,'FaceColor','red','FaceAlpha',.5,'EdgeColor','w','Normalization','probability');
% ax=gca; ax.ThetaTick=[0,90,180,270]; ax.FontWeight='bold'; ax.FontSize=20; ax.RAxisLocation=45;ax.GridAlpha=.5;ax.GridColor='k';
% ax.RTick=[round(linspace(0,.16,4),2)];

% corr(HDc.mean_vector_length(:,:,1),eccentricity.ec(:,:,1))


% [ AllStats ] = ScatterBox(correlationc(:,1,1),correlationt(:,1,1),{'Control','Tilted'},{'SpeedScore'},2)

% [done]=plotstat(controlresult,tiltedresult);


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
%     if fileID==-1
%         disp(['Failed to open',spikefile{i}])
%         continue;
%     end
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
        [SmoothRateMap,nBinsx,nBinsy,occ,Coherence,coherenceJenni]=bindataYoder(Session,60,Session(Session(:,5)==1,:),'no',61);
                
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
        
        results(ii,:,j)=[PeakRate,nSpikes,OverallFR,NumbActiveBins,sparsity,InformationContent,Coherence,Field2Wall,borderScore,FieldWidth,infieldFR,outfieldFR,E,c,p,coherenceJenni];
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


%% PLOT
function done=plotstat(Group1,Group2)
VarNames={'Velocity'};
for c=1
    shadded_line_fig=figure; shadded_line_fig.Color=[1 1 1];shadded_line_fig.OuterPosition=[680 630 729 421];
    % CONTROL
    y=[nanmean(Group1(:,c,1)),nanmean(Group1(:,c,2)),nanmean(Group1(:,c,3)),nanmean(Group1(:,c,4)),nanmean(Group1(:,c,5))];
    SEM=[nanstd(Group1(:,c,1))/sqrt(size(Group1(:,c,1),1)),nanstd(Group1(:,c,2))/sqrt(size(Group1(:,c,2),1)),...
        nanstd(Group1(:,c,3))/sqrt(size(Group1(:,c,3),1)),nanstd(Group1(:,c,4))/sqrt(size(Group1(:,c,4),1)),...
        nanstd(Group1(:,c,5))/sqrt(size(Group1(:,c,5),1))];
    h1=shadedErrorBar(1:5,y,SEM,'-k',0);
    ylabel(VarNames(c));
    ax1=gca; set(ax1,'XTick',[1 2 3 4 5],'Box','off','FontSize',20,'FontWeight','bold','LineWidth',2)
    
    hold on
    
    % TILTED
    y=[nanmean(Group2(:,c,1)),nanmean(Group2(:,c,2)),nanmean(Group2(:,c,3)),nanmean(Group2(:,c,4)),nanmean(Group2(:,c,5))];
    SEM=[nanstd(Group2(:,c,1))/sqrt(size(Group2(:,c,1),1)),nanstd(Group2(:,c,2))/sqrt(size(Group2(:,c,2),1)),...
        nanstd(Group2(:,c,3))/sqrt(size(Group2(:,c,3),1)),nanstd(Group2(:,c,4))/sqrt(size(Group2(:,c,4),1)),...
        nanstd(Group2(:,c,5))/sqrt(size(Group2(:,c,5),1))];
    h2=shadedErrorBar(1:5,y,SEM,'-r',1);
    set(h2.mainLine,'Color',[1 0 0])
    
    %     print(shadded_line_fig, '-dpdf', '-r300',{strcat(FigureLocation,filesep,VarNames(c),'_shadded_line_fig')})
    %     print(strcat(FigureLocation,filesep,VarNames(c),'_shadded_line_fig.pdf'),'-dpdf','-r300','-bestfit')
    %     print(shadded_line_fig,'-bestfit', '-dpdf', '-r600',char(strcat(FigureLocation,filesep,VarNames(c),'_shadded_line_fig')))
    
    %     close all
end

% STATS
for i=1:length(VarNames)
    groups=data.textdata.Session1; groups(1,:)=[]; groups=groups(~cellfun(@isempty,groups));
    measureVal=[Group1(:,i,1),Group1(:,i,2),Group1(:,i,3),Group1(:,i,4),Group1(:,i,5);Group2(:,i,1),Group2(:,i,2),Group2(:,i,3),Group2(:,i,4),Group2(:,i,5)];
    t=table(groups,measureVal(:,1),measureVal(:,2),measureVal(:,3),measureVal(:,4),measureVal(:,5),'VariableNames',{'species','meas1','meas2','meas3','meas4','meas5'});
    Meas = table([1 2 3 4 5]','VariableNames',{'Measurements'});
    
    % CREATE MODELS
    % OVER ALL SESSIONS
    rm1 = fitrm(t,'meas1-meas5~species','WithinDesign',Meas);
    % SESSIONS 1 2 3 - rotation
    rm2 = fitrm(t,'meas1-meas3~species','WithinDesign',Meas(1:3,1));
    % SESSIONS 1 4 5 - dark
    rm3 = fitrm(t,'meas1,meas4,meas5~species','WithinDesign',Meas([1,4,5],1));
    % SESSIONS 1 3 5 - stability of standard sessions
    rm4 = fitrm(t,'meas1,meas3,meas5~species','WithinDesign',Meas([1,3,5],1));
    
    % tbl = mauchly(rm)
    % ranovatbl = ranova(rm)
    disp(strcat('--------------------------------------------------',VarNames(i),'--------------------------------------------------'))
    disp('Multivariate Tests *ALL SESSIONS*')
    manovatbl=manova(rm1)
    disp('Tests of Between-Subjects Effects *ALL SESSIONS*')
    anovatbl=anova(rm1)
    
    disp('Multivariate Tests *ROTATION*')
    manovatbl=manova(rm2)
    disp('Tests of Between-Subjects Effects *ROTATION*')
    anova(rm2)
    
    disp('Multivariate Tests *DARK*')
    manovatbl=manova(rm3)
    disp('Tests of Between-Subjects Effects *DARK*')
    anova(rm3)
    
    disp('Multivariate Tests *STABILITY*')
    manovatbl=manova(rm4)
    disp('Tests of Between-Subjects Effects *STABILITY*')
    anova(rm4)
end

% MEAN AND SEM
for v=1:length(VarNames)
    disp(VarNames(v))
    disp(['Control  S1: ',num2str(nanmean(Group1(:,v,1))),' ± ',num2str(nanstd(Group1(:,v,1)/sqrt(length(Group1(:,v,1))))),...
        ', S2: ',num2str(nanmean(Group1(:,v,2))),' ± ',num2str(nanstd(Group1(:,v,2)/sqrt(length(Group1(:,v,2))))),...
        ', S3: ',num2str(nanmean(Group1(:,v,3))),' ± ',num2str(nanstd(Group1(:,v,3)/sqrt(length(Group1(:,v,3))))),...
        ', S4: ',num2str(nanmean(Group1(:,v,4))),' ± ',num2str(nanstd(Group1(:,v,4)/sqrt(length(Group1(:,v,4))))),...
        ', S5: ',num2str(nanmean(Group1(:,v,5))),' ± ',num2str(nanstd(Group1(:,v,5)/sqrt(length(Group1(:,v,5)))))])
    
    disp(['PAE  S1: ',num2str(nanmean(Group2(:,v,1))),' ± ',num2str(nanstd(Group2(:,v,1)/sqrt(length(Group2(:,v,1))))),...
        ', S2: ',num2str(nanmean(Group2(:,v,2))),' ± ',num2str(nanstd(Group2(:,v,2)/sqrt(length(Group2(:,v,2))))),...
        ', S3: ',num2str(nanmean(Group2(:,v,3))),' ± ',num2str(nanstd(Group2(:,v,3)/sqrt(length(Group2(:,v,3))))),...
        ', S4: ',num2str(nanmean(Group2(:,v,4))),' ± ',num2str(nanstd(Group2(:,v,4)/sqrt(length(Group2(:,v,4))))),...
        ', S5: ',num2str(nanmean(Group2(:,v,5))),' ± ',num2str(nanstd(Group2(:,v,5)/sqrt(length(Group2(:,v,5)))))])
end
done=1;
end


%         try
%             [c,p]=corr(Session(Session(:,5)==1,4),Session(Session(:,5)==1,6));
%             [mean_vector_length,peak_Firing_Rate,preferred_Direction,~,~,Direct_infoContent,~,~,~] = HDCell(Session(:,5),Session(:,7),60);
%         catch
%             c=NaN;
%             p=NaN;
%             mean_vector_length=NaN;
%             peak_Firing_Rate=NaN;
%             preferred_Direction=NaN;
%             Direct_infoContent=NaN;
%         end
%         HD.mean_vector_length(ii,:,j)=mean_vector_length;
%         HD.peak_Firing_Rate(ii,:,j)=peak_Firing_Rate;
%         HD.preferred_Direction(ii,:,j)=preferred_Direction(1);
%         HD.Direct_infoContent(ii,:,j)=Direct_infoContent;


%         meanvel(ii,:,j)=mean(TSXY(Start(j):End(j),4));


%         % bring in ratemap to calc center of field
%         Sess=paths(j:5:length(paths));
%         split=strsplit(Sess{i},'.');
%         matfiles=strcat(split(1),'.',split(2),'_Data.mat');
%         if exist(matfiles{1},'file')~=2;continue;end
%         load(matfiles{1},'field')
%         [row,col]=find(field);
%         k=boundary(row,col);
%         field_XY=[col(k),row(k)];
%         xmax=max(col);xmin=min(col);
%         ymax=max(row);ymin=min(row);
%         x=round(median(xmin:xmax));
%         y=round(median(ymin:ymax));
%
%         pixel=61/max(size(field));
%         x=x*pixel;
%         y=y*pixel;
%
%         xmin=min(Session(:,2));
%         ymax=max(Session(:,3));
%         ymax-y
%         xmin+x
%
%         x=repmat(x,length(Session(:,2)),1);
%         y=repmat(y,length(Session(:,2)),1);
%
%         D2F=sqrt((x-Session(:,2)).^2+(y-Session(:,3)).^2);



%             spkvel=Session(Session(:,5)==1,4);
%             spkifr=Session(Session(:,5)==1,6);
% %             spkd2f=D2F(Session(:,5)==1,:);
%
%             if kstest(spkvel)==1
%                 spkvel=log10(spkvel);
%             end
%             if kstest(spkifr)==1
%                 spkifr=log10(spkifr);
%             end
%             if kstest(spkd2f)==1
%                 spkd2f=log10(spkd2f);
%             end
%
%             tbl=table(spkvel,spkifr,spkifr,'VariableNames',{'VEL','Distance2Field','IFR'});
%             fitlm(tbl,'IFR~Distance2Field+VEL');
%             mdl=fitlm(spkvel,spkifr,'quadratic');
%             rsquared=mdl.Rsquared.Ordinary;

% Tiltedresults=[];
% for i=1:length(Sess1tilt)
%     try
%         % Import video file
%         split=strsplit(Sess1tilt{i},'ASCII');
%         disp(split(1))
%         fileID=fopen(char(strcat(split(1),'VT1_DupRem_PrcFields_lights.txt')),'r');
%         Intro = textscan(fileID,'%s%s%s');
%         fclose(fileID);
%
%         X = [Intro{2}];
%         Y = [Intro{3}];
%
%         X=X(~cellfun('isempty',X));
%         Y=Y(~cellfun('isempty',Y));
%
%         VIDEO=regexprep([X(4:end,1),Y(3:end,1)],',','','emptymatch');
%         N = sscanf(CStr2String(VIDEO, '*'), '%f*');
%         B = reshape(N,[],2);
%
%
%         % interp and smooth
%         %         video(isnan(video))=0;
%         %         [x,y] = InterpolateTrackerZeros_v2(video(:,2),video(:,3));
%         smoothVideo=[SmoothVector(B(:,1),3),SmoothVector(B(:,2),3)];
%
%         % calculate velo
%         [vel_cmPerSec,vel_abs,pixelDist] = InstaVel(smoothVideo,'no',61,60);
%         vel_cmPerSec(vel_cmPerSec>70)=[]; vel_abs(vel_abs>70/pixelDist/60)=[];
%
%         ExtractedAngle=XYangle(smoothVideo(:,1),smoothVideo(:,2));
%
%         % EXTRACT BASIC MOVEMENT DATA
%         AverageAnglePerSec=rad2deg(circ_mean((abs(diff(deg2rad(ExtractedAngle))))*60));
%         OverallDistTraveled=sum(vel_abs*pixelDist);
%         MeanVelocity=mean(vel_cmPerSec);
%         STDVelocity=std(vel_cmPerSec);
%         medVelocity=median(vel_cmPerSec);
%
%         Tiltedresults=[Tiltedresults;AverageAnglePerSec,OverallDistTraveled,MeanVelocity,STDVelocity,medVelocity];
%     catch
%         disp('ERROR SOMEWHERE...DARN IT''S OKAY')
%     end
% end



%
%             for t=1:length(Session)
%                figure(1);plot(Session(t:t+1,2),Session(t:t+1,3));hold on;pause(.0001)
%                if Session(t,4)==1
%                     figure(1);scatter(Session(t,2),Session(t,3),'*r');pause(.0001)
%                end
%
%                figure(2);scatter(t,Session(t,5));title('FR');hold on;pause(.0001)
%
%                figure(3);scatter(t,Session(t,6));title('Vel');hold on;pause(.0001)
%             end

% VarNames=regexprep(VarNames,'R','','emptymatch');
%
% fid = fopen(filename,'r');
% MyTextFile = textscan(fid,'%s','delimiter',',');
% fclose(fid);
% MyTextFile = [MyTextFile{:}];
%
% MyTextFile=MyTextFile(9:end,1);
% MyTextFile(strcmp(MyTextFile,'R'))=[];
%
% eval(MyTextFile{2});
% eval(MyTextFile{3});
% % extract the part containg the matrix to read
% MatrixLines = MyTextFile(9:end);
% % get rid of Chirp1, Chirp2, ...
% MatrixLines = regexprep(MatrixLines,'Chirp[^%d]','');
% % Convet the char into a numerical matrix.
% MyMatrix = str2num(cell2mat(MatrixLines));
% %
%
% fid   = fopen(filename);
% line2 = textscan(fid, '%f%f%f%f%f%f%f%f%f\r\n %*[^\n]','HeaderLines',6);
% fclose(fid);



% MyTextFile=mat(MyTextFile);
%
% MyTextFile{MyTextFile>1000}=[];


% yyaxis left
% scatter(Session(Session(:,5)==1,6), Session(Session(:,5)==1,4))
% yyaxis right
% scatter(Session(Session(:,5)==1,6), D2F(Session(:,5)==1,:))
%
%
%


