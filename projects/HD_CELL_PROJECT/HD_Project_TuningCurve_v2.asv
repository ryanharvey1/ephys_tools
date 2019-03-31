% HD_Project_TuningCurve_v2
% opens .r text files containing HD cell data and creates tuning curves
%
%
% Ryan E Harvey 2018
%
% cd to data and get file names
clear;clc;close all
com=which('HD_Project_TuningCurve_v2');
com=strsplit(com,filesep);
basedir=[com{1},filesep,'Users',filesep,com{3},filesep,'GoogleDrive',filesep,'MatlabDir'];
addpath([basedir,filesep,filesep,'CircStat2012a'],...
    [basedir,filesep,'BClarkToolbox',filesep, 'Analysis']);

if ismac
    path='/Users/ryanharvey/Downloads/HeadDIrectionCells_LauraRyan';
else
    path='D:\Projects\Multi_Region_HD\HeadDIrectionCells_LauraRyan';
end
% path='d:\Users\BClarkLab\Google Drive (lberkowitz@unm.edu)\Manuscripts\In Progress\Ben_HDCProject\Data';
cd(path)
files=dir( '**/*.r');
folders=unique({files.folder});
folders=folders(~ismember(folders,{'Timestamps','timestamps','__MACOSX'}));
folders=folders(~contains(folders,["__MACOSX","Timestamps"]));

for i=1:length(folders)
    area=strsplit(folders{i},filesep);
    area=area{end};
    area=erase(area,'_SpikePositionFiles');
    areas{i}=erase(area,'_');
end
clear area files path

for a=1:length(folders)
    cd(folders{a})
    filenames=dir('*.r');
    filenames={filenames.name}';
    filenames(contains(filenames,'._'))=[];
    
    for i=1:length(filenames)
        % OPEN TEXT FILE AND CLEAN DATA
        disp(['Running:  ',folders{a},filenames{i}])
        fileID=fopen(filenames{i},'r');
        dataArray=textscan(fileID,'%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]','Delimiter',...
            '\t','TextType','string','EmptyValue',NaN,'HeaderLines',2-1,'ReturnOnError',false,'EndOfLine','\r\n');
        fclose(fileID);
        frames=[dataArray{1:end-1}];
        clear dataArray fileID
        % remove last column (all NaNs)
        frames(:,end)=[];
        
        if sum(isnan(frames(:,2)))==size(frames,1)
            continue
        end
        
        HDdata.(areas{a}).id{i,1}=filenames{i};
        
        % convert ts to secs
        samplerate=60;
        frames(:,1)=linspace(0,(length(frames(:,1))/samplerate),length(frames(:,1)));
        
        truesamplerateTS=frames(:,1);
        
        % remove non-detects in the first 4 xy columns
        frames(frames(:,2)==255 | frames(:,3)==255 | frames(:,4)==255 | frames(:,5)==255,:)=[];
        
        % expand frames to create spike binary
        spksover60hz=frames(frames(:,6)>1,:);
        
        idx=find(frames(:,6)>1);
        ts=[];
        for ii=1:size(spksover60hz,1)
            if idx(ii)+1>length(frames)
                tstemp=linspace(frames(idx(ii)-1,1),frames(idx(ii),1),frames(idx(ii),6)+2);
            elseif idx(ii)-1<length(frames)
                tstemp=linspace(frames(idx(ii),1),frames(idx(ii)+1,1),frames(idx(ii),6)+2);
            else
                tstemp=linspace(frames(idx(ii)-1,1),frames(idx(ii)+1,1),frames(idx(ii),6)+2);
            end
            tstemp(1)=[];
            tstemp(end)=[];
            ts=[ts;tstemp'];
        end
        clear tstemp idx
        ts=sort([frames(frames(:,6)<=1,1);ts]);
        
        % remove frames with > spike binary 1 for later recombining
        tempframes=frames;
        tempframes(tempframes(:,6)>1,:)=[];
        
        EXP=zeros(1,size(spksover60hz,2));
        for ii=1:size(spksover60hz,1)
            EXP=[EXP;repmat(spksover60hz(ii,:),spksover60hz(ii,6),1)];
        end
        EXP(1,:)=[];
        EXP(:,6)=ones(size(EXP,1),1);
        
        framesEXP=[tempframes;EXP];
        
        [~,I]=sort(framesEXP(:,1));
        
        framesEXP=framesEXP(I,:);
        
        framesEXP(:,1)=ts;
        
        HDdata.(areas{a}).frames{i,1}=framesEXP;
        
        clear tempframes EXP spksover60hz ii I
        
        %% ////////////////////// SPATIAL ANALYSIS ////////////////////////////
        %
        % BIN SPATIAL RATEMAP
        if contains(areas{a},'ATN')
            arenasize=76;
        else
            arenasize=120;
        end
        binsize=5;
        
        % find center of rat's head
        x=median([framesEXP(:,2),framesEXP(:,4)],2);
        y=median([framesEXP(:,3),framesEXP(:,5)],2);
        
        xedges=linspace(min(x),max(x),arenasize/binsize);
        yedges=linspace(min(y),max(y),arenasize/binsize);
        
        occ=histcounts2(x(framesEXP(:,6)==0),y(framesEXP(:,6)==0),xedges,yedges)/samplerate;
%         occ(occ<.01)=0;
        
        spk=histcounts2(x(framesEXP(:,6)==1),y(framesEXP(:,6)==1),xedges,yedges);
        
        ratemap=(spk./occ);
        
        ratemap(isinf(ratemap))=0;
        
        filtWidth = [3 3]; filtSigma = 1;
        imageFilter=fspecial('gaussian',filtWidth,filtSigma);
        SmoothRateMap = nanconv(ratemap,imageFilter, 'nanout');
        
        SmoothRateMap=rot90(fliplr(SmoothRateMap));
                
        HDdata.(areas{a}).ratemap{i,1}=SmoothRateMap;
        
        [BW,maskedImage,x,y,fieldarea,X] = segmentImage('map',inpaint_nans(SmoothRateMap,2),'binsize',binsize);
        HDdata.(areas{a}).nfields(i,1)=length(fieldarea);
        HDdata.(areas{a}).fieldarea{i,1}=fieldarea;


        
        
        %                     fig=figure;fig.Color=[1 1 1];
        %                     imAlpha=ones(size(SmoothRateMap));
        %                     imAlpha(isnan(SmoothRateMap))=0;
        %                     subplot(1,2,1)
        %                     imagesc(SmoothRateMap,'AlphaData',imAlpha);
        %                     axis xy; colormap jet; axis off; hold on; box off; axis image;
        %
        %                     subplot(1,2,2)
        %                     plot(x(framesEXP(:,6)==0),y(framesEXP(:,6)==0),'.k');hold on
        %                     scatter(x(framesEXP(:,6)==1),y(framesEXP(:,6)==1),'r')
        %                     axis image
        %                     box off
        
        clear spk ratemap filtWidth filtSigma imageFilter binsize
        %
        %% BOARDER ANALYSIS
        % border score
%         [HDdata.(areas{a}).borderscore(i,1),HDdata.(areas{a}).bins2wall(i,:)]=BorderScore(SmoothRateMap);
%         
%         % Border Modulation [Peyrache, Schieferstein, Buzsaki, 2017, DOI:10.1038/s41467-017-01908-3]
%         HDdata.(areas{a}).bordermodulation(i,:)=bordermodulation(length(SmoothRateMap),...
%             framesEXP(:,1),rescale(x,1,length(SmoothRateMap)),rescale(y,1,length(SmoothRateMap)),samplerate,logical(framesEXP(:,6)),3);
%         
%        HDdata.(areas{a}).egocentricmodulation(i,:)=egocentricmodulation(x,...
%            y,framesEXP(:,10),framesEXP,length(SmoothRateMap),samplerate,3);
%         
%         clear top_cor_coef right_cor_coef bottom_cor_coef left_cor_coef...
%             b bywall ifr nframes N edges videots spikets nframes topin...
%             rightin bottomin leftin
%         
        %% GRID SCORE
        gridout=GridScore_Sinusoidal(SpatialAutoCorr(SmoothRateMap,length(SmoothRateMap)),length(SmoothRateMap));
        HDdata.(areas{a}).gridscore(i,1)=gridout.maxSinuGrid;
        
        %% INFO CONTENT
        HDdata.(areas{a}).informationContent(i,1)=infocontent(SmoothRateMap,occ);
        
%         %% SHUFFLE SPATIAL SCORES
%         
%         %Create structured cell array so that shuffling can occur at
%         %correct frame rate (accounts for interpolated timestamps that
%         %alters frame rate of data aka data_video_spk.
%         
%                 secIdx=truesamplerateTS(samplerate+1:samplerate:length(frames),1);
%                 [I,row]=ismember(framesEXP(:,1),secIdx);
%                 I=find(I);
%                 I=[1;I;length(framesEXP)];
%                 for s=1:length(I)
%                     if s==length(I)
%                         datacell{s,1}=framesEXP(I(s):I(end),:);
%                         break
%                     end
%                     datacell{s,1}=framesEXP(I(s):I(s+1)-1,:);
%                 end
%         
%                 x=median([framesEXP(:,2),framesEXP(:,4)],2);
%                 y=median([framesEXP(:,3),framesEXP(:,5)],2);
%         
%                 x2=rescale(x,1,length(SmoothRateMap));
%                 y2=rescale(y,1,length(SmoothRateMap));
%         
%                 shuff_max=(framesEXP(end,1)-20);
%         
%                 bispk=[];
%                 for ishuff=1:500
%                     tempcell=circshift(datacell,randi([20 round(shuff_max)]));
%                     for open=1:length(tempcell)
%                         bispk=[bispk;tempcell{open}(:,6)];
%                     end
%                     bispk=logical(bispk);
%         
%                     % create new ratemap
%                     spk=histcounts2(x(bispk),y(bispk),xedges,yedges);
%                     ratemap=(spk./occ);
%                     ratemap(isinf(ratemap))=0;
%                     ratemap = nanconv(ratemap,fspecial('gaussian',[3 3],1), 'nanout');
%                     ratemap=rot90(fliplr(ratemap));
%         
%                     ratemapstore{ishuff,1}=ratemap;
%                     bispkstore{ishuff,1}=bispk;
%                     bispk=[];
%                 end
%         
%                 ts=framesEXP(:,1);
%                 matsize=length(ratemap);
%                 infocontent_shuff=zeros(500,1);
%                 gridscore_shuff=zeros(500,1);
% %                 borderscore_shuff=zeros(500,1);
% %                 bordermodulation_shuff=zeros(500,1);
%         
%                 try % try in parallel first
%                     parfor ishuff=1:length(ratemapstore)
%                         % create info content distribution
%                         infocontent_shuff(ishuff,1)=infocontent(ratemapstore{ishuff},occ);
%         
%                         % create grid score distribution
%                         try
%                             [gridout] = GridScore_Sinusoidal(SpatialAutoCorr(ratemapstore{ishuff},matsize),matsize);
%                             gridscore_shuff(ishuff,1)=gridout.maxSinuGrid;
%                         catch
%                             gridscore_shuff(ishuff,1)=NaN;
%                         end
%                         
%                         % create border score distribution
% %                         [borderscore_shuff(ishuff,1),~]=BorderScore(ratemapstore{ishuff});
%         
%                         % create border modulation distribution
% %                         bordermodulation_shuff(ishuff,1)=max(bordermodulation(ratemapstore{ishuff},ts,x2,y2,samplerate,bispkstore{ishuff}));
%                     end
%                 catch % out of memory can occur, so we can defalt to the slower for loop
%                     for ishuff=1:length(ratemapstore)
%                         % create info content distribution
%                         infocontent_shuff(ishuff,1)=infocontent(ratemapstore{ishuff},occ);
%         
%                         % create grid score distribution
%                         try
%                             [gridout] = GridScore_Sinusoidal(SpatialAutoCorr(ratemapstore{ishuff},matsize),matsize);
%                             gridscore_shuff(ishuff,1)=gridout.maxSinuGrid;
%                         catch
%                             gridscore_shuff(ishuff,1)=NaN;
%                         end
% 
% %                         [gridscore_shuff(ishuff,1),~,~,~]=gridness_rh(SpatialAutoCorr(ratemapstore{ishuff},matsize));
%         
%                         % create border score distribution
% %                         [borderscore_shuff(ishuff,1),~]=BorderScore(ratemapstore{ishuff});
%         
%                         % create border modulation distribution
% %                         bordermodulation_shuff(ishuff,1)=max(bordermodulation(ratemapstore{ishuff},ts,x2,y2,samplerate,bispkstore{ishuff}));
%                     end
%                 end
%                 % Store results from shuffle
%                 HDdata.(areas{a}).infocontent_shuff_decision(i,1)=HDdata.(areas{a}).informationContent(i,1)>=prctile(infocontent_shuff,95);
%                 HDdata.(areas{a}).gridscore_shuff_decision(i,1)=HDdata.(areas{a}).gridscore(i,1)>=prctile(gridscore_shuff,95);
%                 
%                 HDdata.(areas{a}).infocontent_shuff_95(i,1)=prctile(infocontent_shuff,95);
%                 HDdata.(areas{a}).gridscore_shuff_95(i,1)=prctile(gridscore_shuff,95);
%                 HDdata.(areas{a}).infocontent_shuff_99(i,1)=prctile(infocontent_shuff,99);
%                 HDdata.(areas{a}).gridscore_shuff_99(i,1)=prctile(gridscore_shuff,99);
% %                 HDdata.(areas{a}).borderscore_shuff(i,1)=HDdata.(areas{a}).borderscore(i,1)>=prctile(borderscore_shuff,95);
% %                 HDdata.(areas{a}).bordermodulation_shuff(i,1)=HDdata.(areas{a}).bordermodulation(i,1)>=prctile(bordermodulation_shuff,95);
%         
%                 clear bordermodulation_shuff borderscore_shuff gridscore_shuff...
%                     infocontent_shuff matsize ts bispk spk ratemap ratemapstore...
%                     tempcell bispkstore secIdx x y x2 y2 xedges yedges tstemp row...
%                     I datacell s shuff_max occ ishuff open
%         
% %                 save(['D:\Projects\Multi_Region_HD\tempshuff\','_',num2str(i)],'-struct','HDdata','-v7.3')

        %% //////////////// TEMPORAL ANALYSIS /////////////////
        %
        %THETA MODULATION
        % calc theta mod with
%         [HDdata.(areas{a}).thetaindex(i,1),HDdata.(areas{a}).thetapeak(i,1),...
%             HDdata.(areas{a}).cor(i,:),HDdata.(areas{a}).lag(i,:)]=thetamodulation(framesEXP(logical(framesEXP(:,6)==1),1));
        
        %             % Pass Shuff?
        %             tempframes=framesEXP;
        %             parfor ishuff=1:500
        %                 [thetaindexshuff(ishuff,1),~,~,~]=thetamodulation(tempframes(logical(tempframes(randperm(length(tempframes)),6)),1));
        %             end
        %
        %             if HDdata.(areas{a}).thetaindex(i,1)>=prctile(thetaindexshuff,95)
        %                 HDdata.(areas{a}).modulated(i,1)=1;
        %             else
        %                 HDdata.(areas{a}).modulated(i,1)=0;
        %             end
        
        %% SPEED SCORE
        %         ifr=IFR(framesEXP,samplerate);
        %% ////////////// HD CELL ANAYSIS ///////////////////
        %
        % STABILITY
%         [HDdata.(areas{a}).within_Coeff(i,1),within]=within_HDstability(framesEXP,framesEXP(framesEXP(:,6)==0,:),60,10,6);
        
        %% CREATE TUNING CURVE
        % 6 degree bins
        da=pi/30;
        angBins=[da/2:da:2*pi-da/2];
        % Occupancy
        histAng=hist(frames(:,10),angBins);
        % Number of spikes per bin
        spkPerAng=hist(framesEXP(framesEXP(:,6)==1,10),angBins);
        % Tuning
        hdTuning=(spkPerAng./histAng)*60;
        % remove nan & inf
        hdTuning(isnan(hdTuning) | isinf(hdTuning))=0;
        
        clear da histAng spkPerAng
        
        %% COMPUTE R LENGTH
        rlength(i,1)=circ_r(angBins',hdTuning',deg2rad(6));
        
        %% CIRCSHIFT SO PEAK IS IN THE MIDDLE
        [M,I]=max(hdTuning);
        middlebin=round(median(1:length(hdTuning)));
        
        hdTuningshift=circshift(hdTuning,(middlebin-I)-1);
        
        HDdata.(areas{a}).rawTuningCurve(i,:)=hdTuning;
        
        clear M I middlebin
        
        %% Normalize & Collect tuning curve
        ALL_hdTuning(i,:)=rescale(hdTuningshift,0,1);
        %% Save Temp File
        %         save(['D:\Projects\HDtempData',filesep,areas{a},'_iteration_',num2str(i),'.mat'],'HDdata','-v7.3')
    end
    %% SORT AND SAVE
    [~,I]=sort(rlength);
    ALL_hdTuningsorted=ALL_hdTuning(I,:);
    
    HDdata.(areas{a}).TuningCurve=ALL_hdTuningsorted;
    HDdata.(areas{a}).RLength=rlength;
    
    clear ALL_hdTuning rlength
end

