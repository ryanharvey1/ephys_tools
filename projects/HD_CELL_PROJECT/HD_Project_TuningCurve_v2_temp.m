% HD_Project_TuningCurve_v2_temp
%
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
    [basedir,filesep,'BClarkToolbox',filesep, 'Analysis'],...
    [basedir,filesep,'CMBHOME']);


path='/Users/ryanharvey/Downloads/HeadDIrectionCells_LauraRyan';
cd(path)
files=dir( '**/*.r');
folders=unique({files.folder});

for i=1:length(folders)
    area=strsplit(folders{i},filesep);
    area=area{end};
    area=erase(area,'_SpikePositionFiles');
    areas{i}=erase(area,'_');
end
clear area files path

if ~exist('HDdata','var')
    for a=2:length(folders)
        cd(folders{a})
        filenames=dir('*.r');
        filenames={filenames.name}';
        filenames(contains(filenames,'._'))=[];
        
        for i=1:length(filenames)
            % OPEN TEXT FILE AND CLEAN DATA
            disp(['Running:  ',filenames{i}])
            fileID=fopen(filenames{i},'r');
            dataArray=textscan(fileID,'%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]','Delimiter',...
                '\t','TextType','string','EmptyValue',NaN,'HeaderLines',2-1,'ReturnOnError',false,'EndOfLine','\r\n');
            fclose(fileID);
            frames=[dataArray{1:end-1}];
            clear dataArray fileID
            % remove last column (all NaNs)
            frames(:,end)=[];
            
            % convert ts to secs
            samplerate=60;
            frames(:,1)=linspace(0,(length(frames(:,1))/samplerate),length(frames(:,1)));
            
            truesamplerateTS=frames(:,1);
            
            % remove non-detects in the first 4 xy columns
            frames(frames(:,2)==0 | frames(:,3)==0 | frames(:,4)==0 | frames(:,5)==0,:)=[];
            
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
            
            arenasize=120;
            binsize=3;
            
            % find center of rat's head
            x=median([framesEXP(:,2),framesEXP(:,4)],2);
            y=median([framesEXP(:,3),framesEXP(:,5)],2);
            
            xedges=linspace(min(x),max(x),arenasize/binsize);
            yedges=linspace(min(y),max(y),arenasize/binsize);
            
            occ=histcounts2(x(framesEXP(:,6)==0),y(framesEXP(:,6)==0),xedges,yedges)/samplerate;
            occ(occ<.1)=0;
            
            spk=histcounts2(x(framesEXP(:,6)==1),y(framesEXP(:,6)==1),xedges,yedges);
            
            ratemap=(spk./occ);
            
            ratemap(isinf(ratemap))=0;
            
            filtWidth = [3 3]; filtSigma = 1;
            imageFilter=fspecial('gaussian',filtWidth,filtSigma);
            SmoothRateMap = nanconv(ratemap,imageFilter, 'nanout');
            
            SmoothRateMap=rot90(fliplr(SmoothRateMap));
            
            HDdata.(areas{a}).ratemap{i,1}=SmoothRateMap;
            
            %             fig=figure;fig.Color=[1 1 1];
            %             imAlpha=ones(size(SmoothRateMap));
            %             imAlpha(isnan(SmoothRateMap))=0;
            %             subplot(1,2,1)
            %             imagesc(SmoothRateMap,'AlphaData',imAlpha);
            %             axis xy; colormap jet; axis off; hold on; box off; axis image;
            %
            %             subplot(1,2,2)
            %             plot(x(framesEXP(:,6)==0),y(framesEXP(:,6)==0),'.k');hold on
            %             scatter(x(framesEXP(:,6)==1),y(framesEXP(:,6)==1),'r')
            %             axis image
            %             box off
            
            clear spk ratemap filtWidth filtSigma imageFilter binsize
            %
            %% BOARDER ANALYSIS
            % border score
            [HDdata.(areas{a}).borderscore(i,1),HDdata.(areas{a}).bins2wall(i,:)]=BorderScore(SmoothRateMap);
            
            % Border Modulation [Peyrache, Schieferstein, Buzsaki, 2017, DOI:10.1038/s41467-017-01908-3]
            HDdata.(areas{a}).bordermodulation(i,:)=bordermodulation(SmoothRateMap,...
                framesEXP(:,1),rescale(x,1,length(SmoothRateMap)),rescale(y,1,length(SmoothRateMap)),samplerate,logical(framesEXP(:,6)));

            
            % Egocentric Modulation [Peyrache, Schieferstein, Buzsaki, 2017, DOI:10.1038/s41467-017-01908-3]
            x=rescale(x,1,length(SmoothRateMap));
            y=rescale(y,1,length(SmoothRateMap));
            siz=length(SmoothRateMap);
                        innerboundary=[[6,siz-5];[siz-5,siz-5];[siz-5,6];[6,6];[6,siz-5]];
                        figure;plot(x,y);hold on
                        plot(innerboundary(:,1),innerboundary(:,2))
            
                        in = inpolygon(x,y,innerboundary(:,1),innerboundary(:,2));
                        xinzone=x(~in);
                        yinzone=y(~in);
            
            
                        angles=rad2deg(framesEXP(~in,10));
            
                        boundary=[[linspace(1,siz,100)',repmat(siz,1,100)'];...
                            [repmat(siz,1,100)',linspace(siz,1,100)'];...
                            [linspace(siz,1,100)',ones(1,100)'];...
                            [ones(1,100)',linspace(1,siz,100)']];
            
                        clear left right
                        for j=1:length(xinzone) % calculates min dist from every active bin to the walls
                            for jj=1:length(boundary)
                                d2(jj) = sqrt((xinzone(j)-boundary(jj,1))^2+(yinzone(j)-boundary(jj,2))^2);
                            end
                            [~,I]=min(d2);
                            closestwall=boundary(I,:);
                            towall=atan2d(closestwall(2)-yinzone(j),closestwall(1)-xinzone(j)) + 360*((closestwall(2)-yinzone(j))<0);
                            difference=angles(j)-towall;
            
                            left(j,1)=difference<0 && difference>=-60;
                            right(j,1)=difference>0 && difference<60;
                        end
                        scatter(xinzone(logical(left)),yinzone(logical(left)),'*y')
                        scatter(xinzone(logical(right)),yinzone(logical(right)),'*r')
            
            
                        
                        
                        
                        angle1 = angles+60;
                        angle2 = angles-60;
            
                        x1 = xinzone + 20 * cosd(angles);
                        y1 = yinzone + 20 * sind(angles);
                        
                        x2 = xinzone + 20 * cosd(angle1);
                        y2 = yinzone + 20 * sind(angle1);
                        
                        x3 = xinzone + 20 * cosd(angle2);
                        y3 = yinzone + 20 * sind(angle2);
            %
                        for frames=1:length(xinzone)
                            scatter(xinzone(frames),yinzone(frames),'*r')
                            
                            field=[[x2(frames),y2(frames)];[xinzone(frames),yinzone(frames)];[x(frames),y(frames)]];
                            
                            plot(field(:,1),field(:,2))
                            pause(.01)
                        end
            %
            %             inGC=inpolygon(X,Y,x',y');
            %             GC=(sum(inGC)/size(X,1));
            
            
            clear top_cor_coef right_cor_coef bottom_cor_coef left_cor_coef...
                b bywall ifr nframes N edges videots spikets nframes topin...
                rightin bottomin leftin
            
            %% GRID SCORE
            [HDdata.(areas{a}).gridscore(i,1),~,~,~]=gridness_rh(SpatialAutoCorr(SmoothRateMap,length(SmoothRateMap)));
            
            %% INFO CONTENT
            HDdata.(areas{a}).informationContent(i,1)=infocontent(SmoothRateMap,occ);
            
            %% SHUFFLE SPATIAL SCORES
            
            %Create structured cell array so that shuffling can occur at
            %correct frame rate (accounts for interpolated timestamps that
            %alters frame rate of data aka data_video_spk.
            secIdx=truesamplerateTS(samplerate+1:samplerate:length(frames),1);
            [I,row]=ismember(framesEXP(:,1),secIdx);
            I=find(I);
            I=[1;I;length(framesEXP)];
            for s=1:length(I)
                if s==length(I)
                    datacell{s,1}=framesEXP(I(s):I(end),:);
                    break
                end
                datacell{s,1}=framesEXP(I(s):I(s+1)-1,:);
            end
            
            x=median([framesEXP(:,2),framesEXP(:,4)],2);
            y=median([framesEXP(:,3),framesEXP(:,5)],2);
            
            x2=rescale(x,1,length(SmoothRateMap));
            y2=rescale(y,1,length(SmoothRateMap));
            
            shuff_max=(framesEXP(end,1)-20);
            
            bispk=[];
            for ishuff=1:500
                tempcell=circshift(datacell,randi([20 round(shuff_max)]));
                for open=1:length(tempcell)
                    bispk=[bispk;tempcell{open}(:,6)];
                end
                bispk=logical(bispk);
                
                % create new ratemap
                spk=histcounts2(x(bispk),y(bispk),xedges,yedges);
                ratemap=(spk./occ);
                ratemap(isinf(ratemap))=0;
                ratemap = nanconv(ratemap,fspecial('gaussian',[3 3],1), 'nanout');
                ratemap=rot90(fliplr(ratemap));
                
                ratemapstore{ishuff,1}=ratemap;
                bispkstore{ishuff,1}=bispk;
                bispk=[];
            end
            
            ts=framesEXP(:,1);
            matsize=length(ratemap);
            infocontent_shuff=zeros(500,1);
            gridscore_shuff=zeros(500,1);
            borderscore_shuff=zeros(500,1);
            bordermodulation_shuff=zeros(500,1);
            parfor ishuff=1:length(ratemapstore)
                % create info content distribution
                infocontent_shuff(ishuff,1)=infocontent(ratemapstore{ishuff},occ);
                
                % create grid score distribution
                [gridscore_shuff(ishuff,1),~,~,~]=gridness_rh(SpatialAutoCorr(ratemapstore{ishuff},matsize));
                
                % create border score distribution
                [borderscore_shuff(ishuff,1),~]=BorderScore(ratemapstore{ishuff});
                
                % create border modulation distribution
                bordermodulation_shuff(ishuff,1)=max(bordermodulation(ratemapstore{ishuff},ts,x2,y2,samplerate,bispkstore{ishuff}));
            end
            % Store results from shuffle
            HDdata.(areas{a}).infocontent_shuff(i,1)=HDdata.(areas{a}).informationContent(i,1)>=prctile(infocontent_shuff,95);
            HDdata.(areas{a}).gridscore_shuff(i,1)=HDdata.(areas{a}).gridscore(i,1)>=prctile(gridscore_shuff,95);
            HDdata.(areas{a}).borderscore_shuff(i,1)=HDdata.(areas{a}).borderscore(i,1)>=prctile(borderscore_shuff,95);
            HDdata.(areas{a}).bordermodulation_shuff(i,1)=HDdata.(areas{a}).bordermodulation(i,1)>=prctile(bordermodulation_shuff,95);
            
            clear bordermodulation_shuff borderscore_shuff gridscore_shuff...
                infocontent_shuff matsize ts bispk spk ratemap ratemapstore...
                tempcell bispkstore secIdx x y x2 y2 xedges yedges tstemp row...
                I datacell s shuff_max occ ishuff open 


            %% //////////////// TEMPORAL ANALYSIS /////////////////
            %
            %THETA MODULATION
            % calc theta mod with
            [HDdata.(areas{a}).thetaindex(i,1),HDdata.(areas{a}).thetapeak(i,1),...
                HDdata.(areas{a}).cor(i,:),HDdata.(areas{a}).lag(i,:)]=thetamodulation(framesEXP(logical(framesEXP(:,6)==1),1));
            
            
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
            [HDdata.(areas{a}).within_Coeff(i,1),within]=within_HDstability(framesEXP,framesEXP(framesEXP(:,6)==0,:),60,10,6);
            
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
            
            
        end
        %% SORT AND SAVE
        [~,I]=sort(rlength);
        ALL_hdTuningsorted=ALL_hdTuning(I,:);
        
        HDdata.(areas{a}).TuningCurve=ALL_hdTuningsorted;
        HDdata.(areas{a}).RLength=rlength;
        
        clear ALL_hdTuning rlength
        
    end
end
%% HD CELL EXAMPLES

close all
fig=figure;fig.Color=[1 1 1];
scatter(HDdata.ATN.RLength,HDdata.ATN.within_Coeff,'filled');hold on
scatter(HDdata.PoSDeep.RLength,HDdata.PoSDeep.within_Coeff,'filled')
scatter(HDdata.PoSSup.RLength,HDdata.PoSSup.within_Coeff,'filled')
scatter(HDdata.MECDeepLayers.RLength,HDdata.MECDeepLayers.within_Coeff,'filled')
scatter(HDdata.MECSupLayers.RLength,HDdata.MECSupLayers.within_Coeff,'filled')
scatter(HDdata.PaSDeep.RLength,HDdata.PaSDeep.within_Coeff,'filled')
scatter(HDdata.PaSSup.RLength,HDdata.PaSSup.within_Coeff,'filled')
set(gca,'box','off','LineWidth',2,'FontWeight','bold','FontSize',20,'TickLength',[0;0])
xlabel('RLength')
ylabel('Stability (r)')
legend(areas,'Location','Best');

% print(fig,'-depsc','-tiff','-loose',['/Users/ryanharvey/Downloads',filesep,'RLengthbyStability_Fig.eps'])

close all
atnR=HDdata.ATN.RLength;
posR=[HDdata.PoSDeep.RLength;HDdata.PoSSup.RLength];
mecR=[HDdata.MECDeepLayers.RLength;HDdata.MECSupLayers.RLength];
pasR=[HDdata.PaSDeep.RLength;HDdata.PaSSup.RLength];

atnS=[HDdata.ATN.within_Coeff];
posS=[HDdata.PoSDeep.within_Coeff;HDdata.PoSSup.within_Coeff];
mecS=[HDdata.MECDeepLayers.within_Coeff;HDdata.MECSupLayers.within_Coeff];
pasS=[HDdata.PaSDeep.within_Coeff;HDdata.PaSSup.within_Coeff];


fig=figure; fig.Color=[1 1 1]; fig.OuterPosition=[1 6 960 1052];
findexamples(atnR,atnS,HDdata.ATN.rawTuningCurve);
print(fig,'-depsc','-tiff','-loose',['/Users/ryanharvey/Downloads',filesep,'atn_Fig.eps'])
close all

fig=figure; fig.Color=[1 1 1]; fig.OuterPosition=[1 6 960 1052];
findexamples(posR,posS,[HDdata.PoSDeep.rawTuningCurve;HDdata.PoSSup.rawTuningCurve]);
print(fig,'-depsc','-tiff','-loose',['/Users/ryanharvey/Downloads',filesep,'pos_Fig.eps'])
close all

fig=figure; fig.Color=[1 1 1]; fig.OuterPosition=[1 6 960 1052];
findexamples(mecR,mecS,[HDdata.MECDeepLayers.rawTuningCurve;HDdata.MECSupLayers.rawTuningCurve]);
print(fig,'-depsc','-tiff','-loose',['/Users/ryanharvey/Downloads',filesep,'mec_Fig.eps'])
close all

fig=figure; fig.Color=[1 1 1]; fig.OuterPosition=[1 6 960 1052];
findexamples(pasS,pasR,[HDdata.PaSDeep.rawTuningCurve;HDdata.PaSSup.rawTuningCurve]);
print(fig,'-depsc','-tiff','-loose',['/Users/ryanharvey/Downloads',filesep,'pas_Fig.eps'])
close all

%% BORDER SCORE
close all
ATN=[HDdata.ATN.borderscore,HDdata.ATN.RLength];
PoS=[[HDdata.PoSDeep.borderscore;HDdata.PoSSup.borderscore],[HDdata.PoSDeep.RLength;HDdata.PoSSup.RLength]];
MEC=[[HDdata.MECDeepLayers.borderscore;HDdata.MECSupLayers.borderscore],[HDdata.MECDeepLayers.RLength;HDdata.MECSupLayers.RLength]];
PaS=[[HDdata.PaSDeep.borderscore;HDdata.PaSSup.borderscore],[HDdata.PaSDeep.RLength;HDdata.PaSSup.RLength]];

fig=figure;fig.Color=[1 1 1];
scatter(ATN(:,1),ATN(:,2),'filled');hold on
scatter(PoS(:,1),PoS(:,2),'filled')
scatter(MEC(:,1),MEC(:,2),'filled')
scatter(PaS(:,1),PaS(:,2),'filled')

set(gca,'box','off','LineWidth',2,'FontWeight','bold','FontSize',20,'TickLength',[0;0])
xlabel('Border Score')
ylabel('RLength')
legend({'ATN','PoS','MEC','PaS'},'Location','Best');
%% BORDER MODULATION

close all;clear data
data.ATN=max(HDdata.ATN.bordermodulation,[],2);
data.PoSDeep=max(HDdata.PoSDeep.bordermodulation,[],2);
data.PoSSup=max(HDdata.PoSSup.bordermodulation,[],2);
data.MECDeepLayers=max(HDdata.MECDeepLayers.bordermodulation,[],2);
data.MECSupLayers=max(HDdata.MECSupLayers.bordermodulation,[],2);
data.PaSDeep=max(HDdata.PaSDeep.bordermodulation,[],2);
data.PaSSup=max(HDdata.PaSSup.bordermodulation,[],2);

fig=figure;fig.Color=[1 1 1];
[freq]=ECDF_plot(data,'Border Modulation');
% print(fig,'-depsc','-tiff','-loose',['/Users/ryanharvey/Downloads',filesep,'bordermodulation.eps'])

%% INFO CONTENT
close all;clear data
data.ATN=HDdata.ATN.informationContent;
data.PoSDeep=HDdata.PoSDeep.informationContent;
data.PoSSup=HDdata.PoSSup.informationContent;
data.MECDeepLayers=HDdata.MECDeepLayers.informationContent;
data.MECSupLayers=HDdata.MECSupLayers.informationContent;
data.PaSDeep=HDdata.PaSDeep.informationContent;
data.PaSSup=HDdata.PaSSup.informationContent;

fig=figure;fig.Color=[1 1 1];
[freq]=ECDF_plot(data,'Information Content (bits/spk)');

% print(fig,'-depsc','-tiff','-loose',['/Users/ryanharvey/Downloads',filesep,'bordermodulation.eps'])
print(fig,'-dmeta',['C:\Users\ryanh\Dropbox\school work\UNM\Lab\Projects\Multi_Region_HDcells',filesep,'InformationContent.emf'])

%% Grid Score
close all;clear data
data.ATN=max(HDdata.ATN.gridscore,[],2);
data.PoSDeep=max(HDdata.PoSDeep.gridscore,[],2);
data.PoSSup=max(HDdata.PoSSup.gridscore,[],2);
data.MECDeepLayers=max(HDdata.MECDeepLayers.gridscore,[],2);
data.MECSupLayers=max(HDdata.MECSupLayers.gridscore,[],2);
data.PaSDeep=max(HDdata.PaSDeep.gridscore,[],2);
data.PaSSup=max(HDdata.PaSSup.gridscore,[],2);

fig=figure;fig.Color=[1 1 1];
[freq]=ECDF_plot(data,'Grid Score');

% print(fig,'-depsc','-tiff','-loose',['/Users/ryanharvey/Downloads',filesep,'bordermodulation.eps'])
print(fig,'-dmeta',['C:\Users\ryanh\Dropbox\school work\UNM\Lab\Projects\Multi_Region_HDcells',filesep,'gridscore.emf'])


%% THETA MODULATION FIGURES
clear data
data.ATN=HDdata.ATN.thetaindex;
data.PoSDeep=HDdata.PoSDeep.thetaindex;
data.PoSSup=HDdata.PoSSup.thetaindex;
data.MECDeepLayers=HDdata.MECDeepLayers.thetaindex;
data.MECSupLayers=HDdata.MECSupLayers.thetaindex;
data.PaSDeep=HDdata.PaSDeep.thetaindex;
data.PaSSup=HDdata.PaSSup.thetaindex;

fig=figure;fig.Color=[1 1 1];
[freq]=ECDF_plot(data,'Theta Index');
% print(fig,'-depsc', '-r600',['/Users/ryanharvey/Downloads',filesep,'ThetaModFiglayer.eps'])

%% THETA MODULATION FIGURES
clear data
data.ATN=HDdata.ATN.thetaindex;
data.PoS=[HDdata.PoSDeep.thetaindex;HDdata.PoSSup.thetaindex];
data.MEC=[HDdata.MECDeepLayers.thetaindex;HDdata.MECSupLayers.thetaindex];
data.PaS=[HDdata.PaSDeep.thetaindex;HDdata.PaSSup.thetaindex];

fig=figure;fig.Color=[1 1 1];
[freq]=ECDF_plot(data,'ThetaIndex');
print(fig,'-depsc', '-r600',['/Users/ryanharvey/Downloads',filesep,'ThetaModFigregion.eps'])

%% THETA MODULATED PROPORTIONS
% by layer
close all
y=[sum(HDdata.MECSupLayers.modulated==1)/length(HDdata.MECSupLayers.modulated),sum(HDdata.MECSupLayers.modulated==0)/length(HDdata.MECSupLayers.modulated);...
    sum(HDdata.PaSSup.modulated==1)/length(HDdata.PaSSup.modulated),sum(HDdata.PaSSup.modulated==0)/length(HDdata.PaSSup.modulated);
    sum(HDdata.MECDeepLayers.modulated==1)/length(HDdata.MECDeepLayers.modulated),sum(HDdata.MECDeepLayers.modulated==0)/length(HDdata.MECDeepLayers.modulated);...
    sum(HDdata.PaSDeep.modulated==1)/length(HDdata.PaSDeep.modulated),sum(HDdata.PaSDeep.modulated==0)/length(HDdata.PaSDeep.modulated);...
    sum(HDdata.PoSDeep.modulated==1)/length(HDdata.PoSDeep.modulated),sum(HDdata.PoSDeep.modulated==0)/length(HDdata.PoSDeep.modulated);...
    sum(HDdata.PoSSup.modulated==1)/length(HDdata.PoSSup.modulated),sum(HDdata.PoSSup.modulated==0)/length(HDdata.PoSSup.modulated);...
    sum(HDdata.ATN.modulated==1)/length(HDdata.ATN.modulated),sum(HDdata.ATN.modulated==0)/length(HDdata.ATN.modulated)];
categories={'MECSupLayers','PaSSup','MECDeepLayers','PaSDeep','PoSDeep','PoSSup','ATN'};

barfiglayer=figure;barfiglayer.Color=[1 1 1];
b=barh(y,'stacked')
xlabel('Proportion')
ylabel('Layer/Region')
set(b,'LineStyle','none');
set(gca,'yticklabel',categories,'box','off','LineWidth',2,'FontWeight','bold','FontSize',20,'TickLength',[0;0])
b(1, 1).FaceColor   = [0.929,  0.694,  0.125];
b(1, 2).FaceColor   = [0.3,  0.3,  0.3];

% by region
% close all
mec=[HDdata.MECSupLayers.modulated;HDdata.MECDeepLayers.modulated];
pas=[HDdata.PaSSup.modulated;HDdata.PaSDeep.modulated];
pos=[HDdata.PoSSup.modulated;HDdata.PoSDeep.modulated];
atn=[HDdata.ATN.modulated];

y=[sum(mec==1)/length(mec),sum(mec==0)/length(mec);...
    sum(pas==1)/length(pas),sum(pas==0)/length(pas);...
    sum(pos==1)/length(pos),sum(pos==0)/length(pos);...
    sum(atn==1)/length(atn),sum(atn==0)/length(atn)];
categories={'MEC','PaS','PoS','ATN'};

barfigregion=figure;barfigregion.Color=[1 1 1];
b=barh(y,'stacked')
xlabel('Proportion')
ylabel('Layer/Region')
set(b,'LineStyle','none');
set(gca,'yticklabel',categories,'box','off','LineWidth',2,'FontWeight','bold','FontSize',20,'TickLength',[0;0])
b(1, 1).FaceColor   = [0.929,  0.694,  0.125];
b(1, 2).FaceColor   = [0.3,  0.3,  0.3];

% print(barfiglayer,'-depsc','-tiff','-loose',['/Users/ryanharvey/Downloads',filesep,'ThetaModulationLayer_bar_Fig.eps'])

% print(barfigregion,'-depsc','-tiff','-loose',['/Users/ryanharvey/Downloads',filesep,'ThetaModulationRegion_bar_Fig.eps'])


%% POP VECTOR SORTED BY RLENGTH BY LAYERS
folders=fieldnames(HDdata);

for a=1:length(folders)
    tuning_Fig=figure; tuning_Fig.Color=[1 1 1];
    imagesc(HDdata.(folders{a}).TuningCurve);hold on
    %     plot(rescale(sum(HDdata.(folders{a}).TuningCurve),1,size(HDdata.(folders{a}).TuningCurve,1)),'w','LineWidth',3)
    axis xy;
    box off;
    ylabel('Cells')
    set(gca,'XTick',[1,30.7,60],'XTickLabel',[-180 0 180],'FontSize',20,'FontWeight','bold')
    title(folders{a})
    print(tuning_Fig,'-depsc', '-r600',['/Users/ryanharvey/Downloads',filesep,folders{a},'tuning_Fig.eps'])
    close all
end
%% DIAG POP VEC BY LAYERS

atn=HDdata.ATN.rawTuningCurve;
pos=[HDdata.PoSSup.rawTuningCurve;HDdata.PoSDeep.rawTuningCurve];
mec=[HDdata.MECSupLayers.rawTuningCurve;HDdata.MECDeepLayers.rawTuningCurve];
pas=[HDdata.PaSSup.rawTuningCurve;HDdata.PaSDeep.rawTuningCurve];

[atn]=arrangenorm(HDdata.ATN.rawTuningCurve);
[PoSDeep]=arrangenorm(HDdata.PoSDeep.rawTuningCurve);
[PoSSup]=arrangenorm(HDdata.PoSSup.rawTuningCurve);
[MECDeepLayers]=arrangenorm(HDdata.MECDeepLayers.rawTuningCurve);
[MECSupLayers]=arrangenorm(HDdata.MECSupLayers.rawTuningCurve);
[PaSDeep]=arrangenorm(HDdata.PaSDeep.rawTuningCurve);
[PaSSup]=arrangenorm(HDdata.PaSSup.rawTuningCurve);


data={atn,PoSDeep,PoSSup,MECDeepLayers,MECSupLayers,PaSDeep,PaSSup};
regions=fieldnames(HDdata);


for a=1:length(data)
    tuning_Fig=figure; tuning_Fig.Color=[1 1 1];
    imagesc(data{a});hold on
    axis xy;
    box off;
    ylabel('Cells')
    set(gca,'XTick',[1,30.7,60],'XTickLabel',[-180 0 180],'FontSize',20,'FontWeight','bold')
    title(regions{a})
    print(tuning_Fig,'-depsc', '-r600',['/Users/ryanharvey/Downloads',filesep,regions{a},'tuning_Fig.eps'])
    close all
end

%% POOL BRAIN REGION & CREATE POP VEC SORTED BY RLENGTH
atn=HDdata.ATN.rawTuningCurve;
pos=[HDdata.PoSSup.rawTuningCurve;HDdata.PoSDeep.rawTuningCurve];
mec=[HDdata.MECSupLayers.rawTuningCurve;HDdata.MECDeepLayers.rawTuningCurve];
pas=[HDdata.PaSSup.rawTuningCurve;HDdata.PaSDeep.rawTuningCurve];

[~,atnri]=sort(HDdata.ATN.RLength);
[~,posri]=sort([HDdata.PoSSup.RLength;HDdata.PoSDeep.RLength]);
[~,mecri]=sort([HDdata.MECSupLayers.RLength;HDdata.MECDeepLayers.RLength]);
[~,pasri]=sort([HDdata.PaSSup.RLength;HDdata.PaSDeep.RLength]);

atn=atn(atnri,:);
pos=pos(posri,:);
mec=mec(mecri,:);
pas=pas(pasri,:);

[atn]=shiftnorm(atn);
[pos]=shiftnorm(pos);
[mec]=shiftnorm(mec);
[pas]=shiftnorm(pas);

data={atn,pos,mec,pas};
regions={'atn','pos','mec','pas'};

for a=1:length(data)
    tuning_Fig=figure; tuning_Fig.Color=[1 1 1];
    imagesc(data{a});hold on
    axis xy;
    box off;
    ylabel('Cells')
    set(gca,'XTick',[1,30.7,60],'XTickLabel',[-180 0 180],'FontSize',20,'FontWeight','bold')
    title(regions{a})
    print(tuning_Fig,'-depsc', '-tiff','-loose',['/Users/ryanharvey/Downloads',filesep,regions{a},'tuning_Fig.eps'])
    close all
end
%% DIAG POP VECTORS BY REGION
atn=HDdata.ATN.rawTuningCurve;
pos=[HDdata.PoSSup.rawTuningCurve;HDdata.PoSDeep.rawTuningCurve];
mec=[HDdata.MECSupLayers.rawTuningCurve;HDdata.MECDeepLayers.rawTuningCurve];
pas=[HDdata.PaSSup.rawTuningCurve;HDdata.PaSDeep.rawTuningCurve];

[atn]=arrangenorm(atn);
[pos]=arrangenorm(pos);
[mec]=arrangenorm(mec);
[pas]=arrangenorm(pas);

data={atn,pos,mec,pas};
regions={'atn','pos','mec','pas'};

for a=1:length(data)
    tuning_Fig=figure; tuning_Fig.Color=[1 1 1];
    imagesc(data{a});hold on
    axis xy;
    box off;
    ylabel('Cells')
    set(gca,'XTick',[1,30.7,60],'XTickLabel',[-180 0 180],'FontSize',20,'FontWeight','bold')
    title(regions{a})
    print(tuning_Fig,'-depsc', '-tiff','-loose',['/Users/ryanharvey/Downloads',filesep,regions{a},'tuning_Fig.eps'])
    close all
end


%% LOCAL FUNCTIONS

% SHIFT PEAK TO CENTER
function [out]=shiftnorm(in)
[~,I]=max(in,[],2);
middlebin=round(median(1:size(in,2)));
for i=1:size(in,1)
    out(i,:)=rescale(circshift(in(i,:),(middlebin-I(i,:))-1,2),0,1);
end
end

% FOR DIAG POP VECTOR
function [matout]=arrangenorm(mat)
[~,I]=max(mat,[],2);
[~,I2]=sort(I);
matout=mat(I2,:);
for i=1:size(matout,1)
    matout(i,:)=rescale(matout(i,:),0,1);
end
end

function findexamples(atnR,atnS,tuning)
medR=median(atnR);
medS=median(atnS);
stdR=std(atnR);
stdS=std(atnS);

% TOP
topI=find(atnR>medR+stdR & atnS>medS+stdS);
topI=topI(randi([1,length(topI)],2,1));

% MIDDLE
middleI=find([atnR<medR+stdR & atnR>medR-stdR] & [atnS<medS+stdS & atnS>medS-stdS]);
middleI=middleI(randi([1,length(middleI)],2,1));

% BOTTOM
bottomI=find(atnR<medR-stdR & atnS<medS-stdS);
bottomI=bottomI(randi([1,length(bottomI)],2,1));

subplot(3,2,1)
plot(tuning(topI(1),:),'k','LineWidth',2)
ylabel('Firing Rate (hz)')
title(['RLength: ',num2str(atnR(topI(1))),' Stability: ',num2str(atnS(topI(1)))])
set(gca,'XTick',[1,30.7,60],'XTickLabel',[-180 0 180],'FontSize',20,'FontWeight','bold','TickLength',[0;0],'box','off','LineWidth',2)

subplot(3,2,2)
plot(tuning(topI(2),:),'k','LineWidth',2)
ylabel('Firing Rate (hz)')
title(['RLength: ',num2str(atnR(topI(2))),' Stability: ',num2str(atnS(topI(2)))])
set(gca,'XTick',[1,30.7,60],'XTickLabel',[-180 0 180],'FontSize',20,'FontWeight','bold','TickLength',[0;0],'box','off','LineWidth',2)

subplot(3,2,3)
plot(tuning(middleI(1),:),'k','LineWidth',2)
ylabel('Firing Rate (hz)')
title(['RLength: ',num2str(atnR(middleI(1))),' Stability: ',num2str(atnS(middleI(1)))])
set(gca,'XTick',[1,30.7,60],'XTickLabel',[-180 0 180],'FontSize',20,'FontWeight','bold','TickLength',[0;0],'box','off','LineWidth',2)

subplot(3,2,4)
plot(tuning(middleI(2),:),'k','LineWidth',2)
ylabel('Firing Rate (hz)')
title(['RLength: ',num2str(atnR(middleI(2))),' Stability: ',num2str(atnS(middleI(2)))])
set(gca,'XTick',[1,30.7,60],'XTickLabel',[-180 0 180],'FontSize',20,'FontWeight','bold','TickLength',[0;0],'box','off','LineWidth',2)

subplot(3,2,5)
plot(tuning(bottomI(1),:),'k','LineWidth',2)
ylabel('Firing Rate (hz)')
title(['RLength: ',num2str(atnR(bottomI(1))),' Stability: ',num2str(atnS(bottomI(1)))])
set(gca,'XTick',[1,30.7,60],'XTickLabel',[-180 0 180],'FontSize',20,'FontWeight','bold','TickLength',[0;0],'box','off','LineWidth',2)

subplot(3,2,6)
plot(tuning(bottomI(2),:),'k','LineWidth',2)
ylabel('Firing Rate (hz)')
title(['RLength: ',num2str(atnR(bottomI(2))),' Stability: ',num2str(atnS(bottomI(2)))])
set(gca,'XTick',[1,30.7,60],'XTickLabel',[-180 0 180],'FontSize',20,'FontWeight','bold','TickLength',[0;0],'box','off','LineWidth',2)

end

% function [out] = SpatialAutoCorr(v,nBins)
% 
% com=which('HD_Project_TuningCurve_v2');
% com=strsplit(com,filesep);
% basedir=[com{1},filesep,'Users',filesep,com{3},filesep,'GoogleDrive',filesep,'MatlabDir'];
% addpath([basedir,filesep,filesep,'BClarkToolbox',filesep,'Analyses',filesep,'spikeCode'],...
%     [basedir,filesep,'BClarkToolbox',filesep, 'Analysis']);
% 
% v(isnan(v))=0; v(isinf(v))=0;
% [c,numberOfOverlapPixels] = normxcorr2_general(v,v);
% 
% % reshape autocorr and #ofpixels arrays so they are a column vector
% ci = reshape(c,((nBins*2)-1)*((nBins*2)-1),1);
% pixelsi = reshape(numberOfOverlapPixels,((nBins*2)-1)*((nBins*2)-1),1);
% 
% % keep only bins with >20 samples
% Twentyci = find(pixelsi < 20);
% ci(Twentyci,:) = NaN;
% out = reshape(ci,(nBins*2)-1,(nBins*2)-1);
% end

% function out=infocontent(SmoothRateMap,occ)
% [nBinsx,nBinsy]=size(SmoothRateMap);
% rY=reshape(SmoothRateMap,nBinsx*nBinsy,1);
% rY(isnan(rY) | isinf(rY))=0;
% occRSHP=reshape(occ,nBinsx*nBinsy,1);
% occSUM=sum(occRSHP);
% pX=occRSHP./occSUM;
% [nBins,nCells]=size(rY);
% relR=rY./kron(ones(nBins,1),pX'*rY);
% log_relR=log2(relR);
% log_relR(isinf(log_relR))=0;
% out=sum(kron(pX,ones(1,nCells)).*relR.*log_relR);
% end

% function [out]=bordermodulation(SmoothRateMap,ts,x,y,samplerate,spkbin)
% 
% siz=length(SmoothRateMap);
% 
% top=[[6,siz];[siz-5,siz];[siz-5,siz-5];[6,siz-5];[6,siz]];
% right=[[siz-5,siz-5];[siz,siz-5];[siz,6];[siz-5,6];[siz-5,siz-5]];
% bottom=[[6,6];[siz-5,6];[siz-5,1];[6,1];[6,6]];
% left=[[1,siz-5];[6,siz-5];[6,6];[1,6];[1,siz-5]];
% 
% topin = inpolygon(x,y,top(:,1),top(:,2));
% rightin = inpolygon(x,y,right(:,1),right(:,2));
% bottomin = inpolygon(x,y,bottom(:,1),bottom(:,2));
% leftin = inpolygon(x,y,left(:,1),left(:,2));
% 
% % bin over 50ms
% nframes=60*.05;
% spikets=ts(spkbin);
% videots=ts(~spkbin);
% edges=videots(1:nframes:end);
% [N,~] = histcounts(spikets,edges);
% % convert to spikes/sec
% N=N*(samplerate/nframes);
% % smooth over 200ms
% nframes=round(60*.2);
% ifr=smooth(N,nframes);
% 
% bywall=histcounts(ts(topin,1),edges)';
% bywall(bywall>0)=1;
% b=glmfit(ifr,bywall,'poisson','link','log');
% top_cor_coef=b(2);
% 
% bywall=histcounts(ts(rightin,1),edges)';
% bywall(bywall>0)=1;
% b=glmfit(ifr,bywall,'poisson','link','log');
% right_cor_coef=b(2);
% 
% bywall=histcounts(ts(bottomin,1),edges)';
% bywall(bywall>0)=1;
% b=glmfit(ifr,bywall,'poisson','link','log');
% bottom_cor_coef=b(2);
% 
% bywall=histcounts(ts(leftin,1),edges)';
% bywall(bywall>0)=1;
% b=glmfit(ifr,bywall,'poisson','link','log');
% left_cor_coef=b(2);
% 
% out=[top_cor_coef,right_cor_coef,bottom_cor_coef,left_cor_coef];
% end