% HD_Project_CrossCorr
area={'atn','pos_deep','pos_sup','mec_deep','mec_sup','pas_deep','pas_sup'};

if ~exist('data','var')
    clear;clc;close all
    addpath('/Users/RyanHarvey/GoogleDrive/MatlabDir/BClarkToolbox/Analysis',...
        '/Users/ryanharvey/GoogleDrive/MatlabDir/BClarkToolbox/Analysis/LFP',...
        '/Users/ryanharvey/GoogleDrive/MatlabDir/BClarkToolbox/Analysis/Visualize')
    
    path='H:\HD_Cell_11_29';
    % path='/Users/lauraberkowitz/Google Drive/SpikeTimeCrossCorrTest';
    
    debugging_fig=0;
    
    cd(path)
    
    files=dir( '**/*.r');
    folders=unique({files.folder});
    for i=1:length(folders)
        folders{i} = folders{i}(1:find(folders{i} == filesep, 1, 'last'));
    end
    folders=unique(folders);
    
    area={'atn','pos_deep','pos_sup','mec_deep','mec_sup','pas_deep','pas_sup'};
    
    % loop through regions
    for a=1:length(folders)
        cd(folders{a})
        
        sessions=dir;
        sessions={sessions.name}';
        sessions=sessions(~ismember(sessions,{'.','..','.DS_Store','._.DS_Store'}));
        
        
        peak_crosscorr_amp_dat=[];
        angulardiff_dat=[];
        wavelet_power_0ms_dat=[];
        wavelet_power_dat=[];
        zhist_dat=[];
        comparisons_dat=[];
        
        % loop through sessions
        for s=1:length(sessions)
            cd([folders{a},sessions{s}])
            
            
            filenames=dir('*.ts_R');
            filenames={filenames.name}';
            filenames=erase(filenames,'.ts_R');
            exten={'.ts_R','.r'};
            dotrproblem=0;
            
            % deal with odd file extensions
            if isempty(filenames)
                filenames=dir('*.ts.r');
                filenames={filenames.name}';
                filenames=erase(filenames,'.ts.r');
                dotrproblem=1;
                exten={'.ts.r','.r'};
            end
            
            if length(filenames)==1
                continue
            end
            
            % loop through cells
            for i=1:length(filenames)
                % check for files
                if exist([filenames{i},'.ts_R'],'file')==0 || exist([filenames{i},'.r'],'file')==0
                    if dotrproblem==1
                        if exist([filenames{i},'.ts.r'],'file')==0 || exist([filenames{i},'.r'],'file')==0
                            continue
                        end
                    else
                        continue
                    end
                end
                % OPEN TEXT FILE AND CLEAN DATA
                % spike timestamps
                disp(['Running:  ',folders{a},filenames{i}])
                fileID = fopen([filenames{i},exten{1}],'r');
                dataArray = textscan(fileID,'%f%*s%[^\n\r]','Delimiter',' ','MultipleDelimsAsOne',...
                    true,'TextType','string','EmptyValue',NaN,'HeaderLines',2-1,'ReturnOnError',false,'EndOfLine','\r\n');
                fclose(fileID);
                ts=[dataArray{1:end-1}];
                clear dataArray fileID
                % FIND AND REMOVE NEGATIVE NUMBERS
                ts(ts<0)=[];
                % STORE TIME STAMPS IN CELL ARRAY
                spiketimes{i}=ts;
                
                %% position data
                fileID=fopen([filenames{i},exten{2}],'r');
                dataArray=textscan(fileID,'%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]','Delimiter',...
                    '\t','TextType','string','EmptyValue',NaN,'HeaderLines',2-1,'ReturnOnError',false,'EndOfLine','\r\n');
                fclose(fileID);
                frames=[dataArray{1:end-1}];
                clear dataArray fileID
                % remove last column (all NaNs)
                frames(:,end)=[];
                %remove non-detects in the first 4 xy columns
                frames(frames(:,2)==0 | frames(:,3)==0 | frames(:,4)==0 | frames(:,5)==0,:)=[];
                
                % expand frames to create spike binary
                spksover60hz=frames(frames(:,6)>1,:);
                
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
                
                clear tempframes EXP spksover60hz ii I
                
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
                
                % STORE HD TUNING
                tuningcurves(i,:)=hdTuning;
                
            end
            
            if ~exist('spiketimes','var')
                continue
            end
            
            % remove empties
            tuningcurves=tuningcurves(~cellfun('isempty',spiketimes),:);
            spiketimes=spiketimes(~cellfun('isempty',spiketimes));
            
            if length(spiketimes)==1
                continue
            end
            
            [peak_crosscorr_amp,angulardiff,wavelet_power,wavelet_power_0ms,zhist,comparisons]=cross_corr_(spiketimes,tuningcurves,debugging_fig);
            
            clear spiketimes tuningcurves
            
            % store outcomes
            peak_crosscorr_amp_dat=[peak_crosscorr_amp_dat;peak_crosscorr_amp];
            angulardiff_dat=[angulardiff_dat;angulardiff];
            wavelet_power_0ms_dat=[wavelet_power_0ms_dat;wavelet_power_0ms];
            wavelet_power_dat=cat(3,wavelet_power_dat,wavelet_power);
            zhist_dat=[zhist_dat;zhist];
            comparisons_dat=[comparisons;comparisons_dat];
            
            if debugging_fig==1
                figure(2)
                axH = findall(gcf,'type','axes');
                set(axH,'clim',[0 max([axH.CLim])])
                
                figure
                scatter(peak_crosscorr_amp,angulardiff)
                
                figure
                plot(mean(wavelet_power,1));
            end
        end
        data.(area{a}).peak_crosscorr_amp=peak_crosscorr_amp_dat;
        data.(area{a}).angulardiff=angulardiff_dat;
        data.(area{a}).wavelet_power_0ms=wavelet_power_0ms_dat;
        data.(area{a}).wavelet_power=wavelet_power_dat;
        data.(area{a}).zhist=zhist_dat;
        data.(area{a}).comparisons=comparisons_dat;
    end
end

fig1=figure(1);fig1.Color=[1 1 1];
fig2=figure(2);fig2.Color=[1 1 1];
fig3=figure(3);fig3.Color=[1 1 1];

co=[ 0    0.4470    0.7410;
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560;
    0.4660    0.6740    0.1880;
    0.3010    0.7450    0.9330;
    0.6350    0.0780    0.1840];
categories=erase(area,'_');
for i=1:length(area)
    % remove nans
    nanidx=sum(isnan(data.(area{i}).wavelet_power_0ms),2)==size(data.(area{i}).wavelet_power_0ms,2);
    wavelet_power_0ms=data.(area{i}).wavelet_power_0ms(~nanidx,:);
    angulardiff=data.(area{i}).angulardiff(~nanidx,:);
    
    figure(1)
    subplot(2,4,i)
    scatter(data.(area{i}).angulardiff,data.(area{i}).peak_crosscorr_amp,'k','filled');
    title(erase(area{i},'_'))
    xlabel('Angular difference (�)')
    ylabel('Peak amplitude (z)')
    box off
    axis tight
    axis square
    xlim([0 180])
    
    figure(2)
    errBar=nanstd(wavelet_power_0ms)/sqrt(size(wavelet_power_0ms,1));
    y=nanmean(wavelet_power_0ms);
    x=1:size(wavelet_power_0ms,2);
    varargout=shadedErrorBar(x,y,errBar,{'-','markerfacecolor',[co(i,:)]},1);
    xlabel('Frequency (Hz)')
    ylabel('Power (dB)')
    box off
    axis tight
    hold on
    
    figure(3)
    subplot(2,4,i)
    [~,I]=max(wavelet_power_0ms,[],2);
    meanFreq_0ms=I*10;
    scatter(angulardiff,meanFreq_0ms,'k','filled');
    title(erase(area{i},'_'))
    xlabel('Angular difference (�)')
    ylabel('Mean Frequency at 0 ms')
    box off
    axis tight
    axis square
    xlim([0 180])
    ylim([0 500])
end

figure(1)
axH = findall(gcf,'type','axes');
set(axH,'ylim',[min([axH.YLim]) max([axH.YLim])],'LineWidth',2,...
    'FontWeight','bold','FontSize',15,'XTick',linspace(0,180,5),'XTickLabel',linspace(0,180,5))

figure(2);
[h,icons]=legend(categories,'FontSize',12,'Location','best');
legend('boxoff')
set(gca,'box','off','LineWidth',2,'FontWeight','bold','FontSize',15)
set(gca,'XTick',linspace(1,50,5),'XTickLabel',linspace(0,500,5))

figure(3)
axH = findall(gcf,'type','axes');
set(axH,'LineWidth',2,'FontWeight','bold','FontSize',15,'XTick',linspace(0,180,5),'XTickLabel',linspace(0,180,5))

% print(figure(2),'-dpng','-r1200',['C:\Users\ryanh\GoogleDrive',filesep,'ljsdflsjdf.png'])


function [peak_crosscorr_amp,angulardiff,wavelet_power,wavelet_power_0ms,zhistsave,comparisons]=cross_corr_(spikearray,tuningcurves,debugging_fig)
% Plots CCG for all cell pairs
%
%   Input:
%           spikearray: cell array of spike times in seconds
%
%
% Ryan Harvey 2018

peak_crosscorr_amp=[];
angulardiff=[];
wavelet_power=[];
wavelet_power_0ms=[];
comparisons=[];
zhistsave=[];


rc=length(spikearray);
index = reshape(1:rc^2, rc, rc).';
for i=1:rc
    for ii=i+1:rc
        zhist=z_xcorr(spikearray{i},spikearray{ii});
        x=-1000:.5:1000;
        restrict=find(abs(x)==20);
        [phase,pow]=multiphasevec(1:10:500,zhist(restrict(1):restrict(2)),2000,6);
        
        if debugging_fig==1
            figure(1)
            subplot(rc,rc,index(ii,i))
            h=area(x(restrict(1):restrict(2)),zhist(restrict(1):restrict(2)));
            h.FaceColor = [0 0 0];
            box off
            grid off
            title([num2str(i),' vs. ',num2str(ii)])
            pause(.0001)
            
            figure(2)
            subplot(rc,rc,index(ii,i))
            imagesc(pow)
            set(gca,'YDir','normal')
            title([num2str(i),' vs. ',num2str(ii)])
            colormap jet
            
            figure(3)
            subplot(rc,rc,index(ii,i))
            plot(rescale(tuningcurves(i,:),0,1),'k')
            hold on
            plot(rescale(tuningcurves(ii,:),0,1),'r')
            set(gca,'XTick',[1,30.7,60],'XTickLabel',[-180 0 180],'TickLength',[0;0],'box','off')
            title([num2str(i),' vs. ',num2str(ii)])
        end
        
        % outcome measures
        % degree offset
        [~,peakI1]=max(tuningcurves(i,:));
        [~,peakI2]=max(tuningcurves(ii,:));
        normDeg=mod(peakI1*6-peakI2*6,360);
        angulardiff=[angulardiff;min(360-normDeg,normDeg)];
        
        % peak cross corr amp at 0ms (mean of +-.5ms)
        restrict=find(abs(x)==.5);
        peak_crosscorr_amp=[peak_crosscorr_amp;nanmean(zhist(restrict(1):restrict(2)))];
        
        % wavelet_power at 0ms
        wavelet_power_0ms=[wavelet_power_0ms;pow(:,[-20:.5:20]==0)'];
        
        %wavelet_power
        wavelet_power=cat(3,wavelet_power,pow);
        
        % comparison ids
        comparisons=[comparisons;[i,ii]];
        
        % cross corr
        zhistsave=[zhistsave;zhist];
        
    end
end
end


