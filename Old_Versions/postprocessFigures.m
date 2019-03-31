function p=postprocessFigures(data,cellid)
com=which('postprocessFigures');
com=strsplit(com,filesep);

basedir=[com{1},filesep,'Users',filesep,com{3},filesep,'GoogleDrive',filesep,'MatlabDir'];
addpath([basedir,filesep,'BClarkToolbox',filesep,'Analysis',filesep,'Visualize',filesep,'panel'],...
    [basedir,filesep,'BClarkToolbox',filesep,'Analysis'],...
    [basedir,filesep,'BClarkToolbox',filesep,'Analysis',filesep,'Utils'],...
    [basedir,filesep,'CircStat2012a'],...
    [basedir,filesep,'plotSpikeRaster']);
    
% how many cells, how many sessions?
ncells=length(data.Spikes);
nsessions=size(data.events,2);
if any(any(contains(data.mazetypes,'CircularTrack') | contains(data.mazetypes,'LinearTrack')))
    nsessions=size(data.events,2)+1;
end

% set up colormap
theta=0:.01:2*pi;
c=hsv(length(theta));

if any(any(contains(data.mazetypes,'CircularTrack') | contains(data.mazetypes,'LinearTrack')))
    unlinear=data.linear_track.nonlinearFrames;
end

if exist('cellid','var')
    cells=find(contains(data.spikesID.TetrodeNum,cellid{1}) & ismember(data.spikesID.CellNum,cellid{2}))';
    plotraster=0;
    xcorr=0;
else
    cells=1:ncells;
    plotraster=1;
    if ncells>20
        xcorr=0;
    else
        xcorr=1;
    end
end


for i=cells
    tetrode=strsplit(data.spikesID.paths{i},filesep);
    tetrode=tetrode{end};
    trodeID=str2double(extractBetween(tetrode,'TT','.'));
    
    % SET UP FIGURE 
    fig=figure('Name',[data.rat,'  ',data.sessionID,'  ',tetrode,' Cell: ',num2str(data.spikesID.CellNum(i))],'NumberTitle','off');
    fig.Color=[1 1 1];
    p = panel(fig); 
    p.pack(7, nsessions); 
    % set margins
    p.de.margin = 7;

    % and some properties
    p.fontsize = 12;
    p.fontweight='bold';

    for ns=1:nsessions
        p(1, ns).select(); 
        ax=gca; 

        if any(contains(data.mazetypes,'CircularTrack') | contains(data.mazetypes,'LinearTrack'))
            
            if ns==1 % FOR RIGHT RUNNING DIRECTION
                for lap=1:length(data.linear_track.right{1,i}.laps)
                    unlinearR{lap,1}=unlinear(ismember(unlinear(:,1),...
                        data.linear_track.right{1,i}.laps{lap}(:,1)),:);
                end
                spkframes=data.linear_track.right{1,i}.dataspks;
                spkts=spkframes(spkframes(:,6)==1,1);
                unlinear_togetherR=vertcat(unlinearR{:});
                if contains(data.mazetypes{1},'LinearTrack')
                    [unlinear_togetherR(:,2),unlinear_togetherR(:,3)]=...
                        ParameterizeToLinearTrack2(unlinear_togetherR(:,2),unlinear_togetherR(:,3));
                end
                plot(unlinear_togetherR(:,2),unlinear_togetherR(:,3),'.k');hold on
                [ts,idx]=unique(unlinear_togetherR(:,1));
                scatter(interp1(ts,unlinear_togetherR(idx,2),spkts),interp1(ts,unlinear_togetherR(idx,3),spkts),10,...
                    interp1(data.lfp.ts,...
                    data.lfp.theta_phase(trodeID,:),spkts,'linear'),'filled');
                box off;axis image; axis off
                colormap(ax,c)
                title(['Left Trials, nSpikes: ',num2str(length(spkts))]);
                
                % PLOT SPIKES ON LINEARIZED PATH WITH EACH SPIKE COLOR
                % CODED BY ITS THETA PHASE
                p(2, ns).select(); 
                ax=gca;
                ts=data.linear_track.right{1,i}.dataspks(:,1);
                x=data.linear_track.right{1,i}.dataspks(:,2);
                x=rescale(x,1,length(data.ratemap{i,ns}));
                spkbinary=logical(data.linear_track.right{1,i}.dataspks(:,6));
                plot(x,ts,'.k');
                axis tight
                hold on;box off; axis off
                scatter(x(spkbinary),ts(spkbinary),20,...
                    interp1(data.lfp.ts,...
                    data.lfp.theta_phase(trodeID,:),ts(spkbinary)','linear'),'filled');
                colormap(ax,c)
                
                % PLOT EACH LAPS RATEMAP WITH THE OVERALL FIRING RATE
                % SUPERIMPOSED
                p(3, ns).select(); 
                ax = gca; 
                laps=reshape([data.linear_track.right{1,i}.maps{:}],...
                    [],length(data.linear_track.right{1,i}.maps));
                imagesc(flipud(imrotate(laps,90)));hold on;set(gca,'Ydir','Normal')
                ratemap=data.ratemap{i,ns};
                plot(rescale(ratemap,1,size(laps,2)),'LineWidth',2, 'color','w');
                set(gca,'XTickLabel',[]);
                axis tight
                hold on;box off;
                colormap(ax,jet(255))
                title(sprintf('Info Content: %4.2f Peak Rate: %4.2f',...
                    data.measures(i,contains(data.varnames,["InformationContent","PeakRate"]),ns)))

                
                % PHASE BY POSITION 
                p(4, ns).select();
                ax = gca;
                x=rescale(x,0,1);
                phase=interp1(data.lfp.ts,...
                    data.lfp.theta_phase(trodeID,:),ts(spkbinary)','linear');
                 scatter([x(spkbinary);x(spkbinary)],[phase';phase'+2*pi]*180/pi,15,'Filled','k');
                box off;
                xlim([min(x) max(x)])
                set(ax,'YTick',0:360:720,'YTickLabel',[{'0'},{'2\pi'},{'4\pi'}])

                if isfield(data.linear_track.right{i},'fields')
                    for f=1:length(data.linear_track.right{1,i}.fields)
                        rho(f)=data.linear_track.right{1,i}.fields{f}.ThPrecess.circLinCorr;
                        pval(f)=data.linear_track.right{1,i}.fields{f}.ThPrecess.pval;
                    end
                    a3=[rho;pval];
                    title(sprintf(repmat(' rho:%4.2f p:%4.2f,',1,length(rho)),a3(:)'))
                    clear rho pval
                    hold on
                    % plot correlation lines
                    plot_corr_line(x(spkbinary),length(ratemap),data.linear_track.right{1,i})
                end
                
                % PHASE MAP
                p(5, ns).select();
                ax = gca;
                ts=data.linear_track.right{1,i}.dataspks(data.linear_track.right{1,i}.dataspks(:,6)==0,1);
                x=data.linear_track.right{1,i}.dataspks(data.linear_track.right{1,i}.dataspks(:,6)==0,2);
                spkts=data.linear_track.right{1,i}.dataspks(data.linear_track.right{1,i}.dataspks(:,6)==1,1);
                [~,I]=unique(ts);
                ts=ts(I);
                x=x(I);
                Phasemap(ts,x,spkts,data.lfp,length(ratemap),trodeID,data.samplerate)
                colormap(ax,jet(255))
                
                % PLOT AUTOCORRELATION
                p(6, ns).select(); 
                ax = gca; 
                y=data.thetaautocorr{i,ns};
                plot(y,'LineWidth',2, 'color','k');
                axis tight
                hold on;box off; axis off
                title(sprintf('Theta modulation: %4.2f',...
                    data.measures(i,contains(data.varnames,'thetaindex'),ns)))
                
                % PLOT AVERAGE WAVEFORMS [make this compatible with less channels]
                p(7, ns).select(); 
                ax = gca; 
                waves=zeros(4,150);
                waves(1:size(data.avgwave{i},1),:)=data.avgwave{i};
                plot(1:150,waves(1,:),'LineWidth',2, 'color','k');hold on
                plot(151:300,waves(2,:),'LineWidth',2, 'color','k');
                plot(1:150,waves(3,:)+max(abs([waves(1,:),waves(2,:)])),'LineWidth',2, 'color','k');
                plot(151:300,waves(4,:)+max(abs([waves(1,:),waves(2,:)])),'LineWidth',2, 'color','k');
                axis tight
                hold on;box off; axis off
                title(sprintf('ShortISI: %4.2f :IsoDist %4.2f',...
                    data.measures(i,contains(data.varnames,["IsolationDistance","ShortISI"]),ns)))
                
            elseif ns==2 % FOR LEFT RUNNING DIRECTION
               
                for lap=1:length(data.linear_track.left{1,i}.laps)
                    unlinearL{lap,1}=unlinear(ismember(unlinear(:,1),...
                        data.linear_track.left{1,i}.laps{lap}(:,1)),:);
                end
                spkframes=data.linear_track.left{1,i}.dataspks;
                spkts=spkframes(spkframes(:,6)==1,1);
                unlinear_togetherL=vertcat(unlinearL{:});
                if contains(data.mazetypes{1},'LinearTrack')
                    [unlinear_togetherL(:,2),unlinear_togetherL(:,3)]=...
                        ParameterizeToLinearTrack2(unlinear_togetherL(:,2),unlinear_togetherL(:,3));
                end
                plot(unlinear_togetherL(:,2),unlinear_togetherL(:,3),'.k');hold on
                [ts,idx]=unique(unlinear_togetherL(:,1));
                scatter(interp1(ts,unlinear_togetherL(idx,2),spkts),interp1(ts,unlinear_togetherL(idx,3),spkts),10,...
                    interp1(data.lfp.ts,...
                    data.lfp.theta_phase(trodeID,:),spkts,'linear'),'filled');
                box off; axis image;axis off
                colormap(ax,c)
                title(['Right Trials, nSpikes: ',num2str(length(spkts))]);

                
                % PLOT SPIKES ON LINEARIZED PATH WITH EACH SPIKE COLOR
                % CODED BY ITS THETA PHASE
                p(2, ns).select(); 
                ax = gca; 
                ts=data.linear_track.left{1,i}.dataspks(:,1);
                x=data.linear_track.left{1,i}.dataspks(:,2);
                x=rescale(x,1,length(data.ratemap{i,ns}));
                spkbinary=logical(data.linear_track.left{1,i}.dataspks(:,6));
                plot(x,ts,'.k');
                axis tight
                hold on;box off; axis off
                scatter(x(spkbinary),ts(spkbinary),20,...
                    interp1(data.lfp.ts,...
                    data.lfp.theta_phase(trodeID,:),ts(spkbinary)','linear'),'filled');
                colormap(ax,c)
%                 title([tetrode,' Cell: ',num2str(data.spikesID.CellNum(i)),'  nSpikes: ',num2str(sum(spkbinary))]);
                
                % PLOT EACH LAPS RATEMAP WITH THE OVERALL FIRING RATE
                % SUPERIMPOSED
                p(3, ns).select(); 
                ax = gca; 
                laps=reshape([ data.linear_track.left{1,i}.maps{:}],...
                    [],length(data.linear_track.left{1,i}.maps));
                imagesc(flipud(imrotate(laps,90)));hold on;set(gca,'Ydir','Normal')
                ratemap= data.ratemap{i,ns};
                plot(rescale(ratemap,1,size(laps,2)),'LineWidth',2, 'color','w');
                set(gca,'XTickLabel',[]);
                axis tight
                hold on;box off;
                colormap(ax,jet(255))
%                 title([num2str(round(max(max(ratemap)),2)),' Hz']);
                title(sprintf('Info Content: %4.2f Peak Rate: %4.2f',...
                    data.measures(i,contains(data.varnames,["InformationContent","PeakRate"]),ns)))

                % PHASE BY POSITION
                p(4, ns).select();
                ax = gca;
                x=rescale(x,0,1);
                phase=interp1(data.lfp.ts,...
                    data.lfp.theta_phase(trodeID,:),ts(spkbinary)','linear');
                 scatter([x(spkbinary);x(spkbinary)],[phase';phase'+2*pi]*180/pi,15,'Filled','k');
                 axis tight;box off;axis off
                 xlim([min(x) max(x)])
                 
                 %                 title(sprintf('Corr: %4.2f p: %4.2f',...
                 %                     data.measures(i,contains(data.varnames,["PhcircLinCorr","Phpval"]),ns)))
                 if isfield(data.linear_track.left{i},'fields')
                     for f=1:length(data.linear_track.left{1,i}.fields)
                         rho(f)=data.linear_track.left{1,i}.fields{f}.ThPrecess.circLinCorr;
                         pval(f)=data.linear_track.left{1,i}.fields{f}.ThPrecess.pval;
                     end
                     a3=[rho;pval];
                     title(sprintf(repmat(' rho:%4.2f p:%4.2f,',1,length(rho)),a3(:)'))
                     clear rho pval
                     
                     hold on
                     % plot correlation lines
                     plot_corr_line(x(spkbinary),length(ratemap),data.linear_track.left{1,i})
                 end
                 
                 p(5, ns).select();
                 ax = gca;
                 ts=data.linear_track.left{1,i}.dataspks(data.linear_track.left{1,i}.dataspks(:,6)==0,1);
                 x=data.linear_track.left{1,i}.dataspks(data.linear_track.left{1,i}.dataspks(:,6)==0,2);
                 spkts=data.linear_track.left{1,i}.dataspks(data.linear_track.left{1,i}.dataspks(:,6)==1,1);
                 [~,I]=unique(ts);
                 ts=ts(I);
                 x=x(I);
                 Phasemap(ts,x,spkts,data.lfp,length(ratemap),trodeID,data.samplerate)
                 colormap(ax,jet(255))
                
                
                % PLOT AUTOCORRELATION
                p(6, ns).select(); 
                ax = gca; 
                y=data.thetaautocorr{i,ns};
                plot(y,'LineWidth',2, 'color','k');
                axis tight
                hold on;box off; axis off
                title(sprintf('Theta modulation: %4.2f',...
                    data.measures(i,contains(data.varnames,'thetaindex'),ns)))
                
                p(7, ns).select();
                ax=gca;
                cb=imagesc(theta);hold on
                colormap(ax,c)
                set(ax,'XTick',linspace(0,length(theta),2),'XTickLabel',...
                    [{'0'} {'2\pi'}],'TickLength',[0,0],'YTick',[.5 1 1.5],'YTickLabel',[])
                x = 0:1:length(theta);
                y = rescale(gaussmf(x,[median([1,length(theta)/2]) length(theta)/2]),.5,1.5);
                plot(x,y,'w')
                axis square
                xlim([(min(x))*2 max(x)*2])
                
            elseif ns>2 % FOR OPEN FIELD AND BOX
                % PLOT TUNNING CURVE
                p(1, ns).select(); 
                [data_video_spk,~]=createframes_w_spikebinary(data,ns-1,i);
                [~,~,~,~,~,tuning]=tuningcurve(data_video_spk(data_video_spk(:,6)==0,4),...
                    data_video_spk(data_video_spk(:,6)==1,4),data.samplerate);
                angBins=0:6:360;
                bin_centers=movmedian(angBins,2);
                bin_centers(1)=[];
                polarplot(deg2rad(bin_centers),tuning,'k')
                ax=gca;
                set(ax,'RGrid','off','ThetaGrid','off','ThetaTick',[0 90],...
                    'ThetaTickLabels',[0 90],'RTick',max(tuning),'RTickLabel',max(tuning),'LineWidth',3);
                axis tight
                
                % PLOT SPIKES ON PATH WITH EACH SPIKE COLOR
                % CODED BY ITS THETA PHASE
                p(2, ns).select(); 
                ax=gca;
                plot(data_video_spk(:,2),data_video_spk(:,3),'LineWidth',1,'color','k');
                hold on; axis off
                scatter(data_video_spk(data_video_spk(:,6)==1,2),data_video_spk(data_video_spk(:,6)==1,3),10,...
                    interp1(data.lfp.ts,...
                    data.lfp.theta_phase(trodeID,:),data_video_spk(data_video_spk(:,6)==1,1),'linear'),'filled');
                box off; axis image
                colormap(ax,c)
                title(['nSpikes: ',num2str(sum(data_video_spk(:,6)==1))]);

                % PLOT RATEMAP
                p(3, ns).select(); 
                ax = gca; 
                SmoothRateMap=data.ratemap{i,ns};
                imAlpha=ones(size(SmoothRateMap));
                imAlpha(isnan(SmoothRateMap))=0;
                imagesc(SmoothRateMap,'AlphaData',imAlpha);
                axis xy; axis off; hold on; box off; axis image;
                colormap(ax,jet(255))
                %                 title([num2str(round(max(max(SmoothRateMap)),2)),' Hz']);
                title(sprintf('Info Content: %4.2f Peak Rate: %4.2f',...
                    data.measures(i,contains(data.varnames,["InformationContent","PeakRate"]),ns)))

                % PHASE BY POSITION 
                p(4, ns).select();
                ax = gca;
                if size(data.ThPrecess{i,ns},1)>1
                    scatter([data.ThPrecess{i,ns}(:,1);...
                        data.ThPrecess{i,ns}(:,1)],...
                        [data.ThPrecess{i,ns}(:,2);...
                        data.ThPrecess{i,ns}(:,2)+360],10,'Filled','k');
                end
                axis tight;box off;axis off;axis square
                title(sprintf('Corr: %4.2f p: %4.2f',...
                    data.measures(i,contains(data.varnames,["PhcircLinCorr","Phpval"]),ns)))
                hold on
                plot_corr_line(data.ThPrecess{i,ns}(:,1),length(ratemap),...
                    data.measures(i,contains(data.varnames,["slopeCpU","phaseOffset"]),ns))
                
%                 p(5, ns).select();
%                 ax = gca;
%                 ts=data_video_spk(data_video_spk(:,6)==0,1);
%                 x=data_video_spk(data_video_spk(:,6)==0,2);
%                 spkts=data_video_spk(data_video_spk(:,6)==1,1);
%                 [~,I]=unique(ts);
%                 ts=ts(I);
%                 x=x(I);
%                 Phasemap(ts,x,spkts,data.lfp,length(ratemap),trodeID,data.samplerate)
%                 colormap(ax,jet(255))
%                 axis square
                
                % PLOT AUTOCORRELATION
                p(5, ns).select(); 
                ax = gca; 
                y=data.thetaautocorr{i,ns};
                plot(y,'LineWidth',2, 'color','k');
                axis tight
                hold on;box off; axis off
                title(sprintf('Theta modulation: %4.2f',...
                    data.measures(i,contains(data.varnames,'thetaindex'),ns)))
            end
            
        elseif contains(data.mazetypes{ns},'Cylinder') % IF FIRST MAZE IS NOT A TRACK
            
            % PLOT TUNNING CURVE
            p(1, ns).select();
            [data_video_spk,~]=createframes_w_spikebinary(data,ns,i);
            [~,~,~,~,~,tuning]=tuningcurve(data_video_spk(data_video_spk(:,6)==0,4),...
                data_video_spk(data_video_spk(:,6)==1,4),data.samplerate);
            angBins=0:6:360;
            bin_centers=movmedian(angBins,2);
            bin_centers(1)=[];
            polarplot(deg2rad(bin_centers),tuning,'k')
            ax=gca;
            set(ax,'RGrid','off','ThetaGrid','off','ThetaTick',[0 90],...
                'ThetaTickLabels',[0 90],'RTick',max(tuning),'RTickLabel',max(tuning),'LineWidth',3);
            axis tight
            
            % PLOT SPIKES ON LINEARIZED PATH WITH EACH SPIKE COLOR
            % CODED BY ITS THETA PHASE
            p(2, ns).select(); 
            ax=gca;
            [data_video_spk,~]=createframes_w_spikebinary(data,ns,i);
            plot(data_video_spk(:,2),data_video_spk(:,3),'LineWidth',1,'color','k');
            hold on; axis off
            scatter(data_video_spk(data_video_spk(:,6)==1,2),...
                data_video_spk(data_video_spk(:,6)==1,3),10,...
                interp1(data.lfp.ts,...
                data.lfp.theta_phase(trodeID,:),data_video_spk(data_video_spk(:,6)==1,1),'linear'),'filled');
            box off; axis image
            colormap(ax,c)
            title(['nSpikes: ',num2str(sum(data_video_spk(:,6)==1))]);
            
            % PLOT RATEMAP
            p(3, ns).select();
            ax=gca;
            SmoothRateMap=data.ratemap{i,ns};
            imAlpha=ones(size(SmoothRateMap));
            imAlpha(isnan(SmoothRateMap))=0;
            imagesc(SmoothRateMap,'AlphaData',imAlpha);
            axis xy; axis off; hold on; box off; axis image;
            colormap(ax,jet(255))
            %             title([num2str(round(max(max(SmoothRateMap)),2)),' Hz']);
            title(sprintf('Info Content: %4.2f Peak Rate: %4.2f',...
                data.measures(i,contains(data.varnames,["InformationContent","PeakRate"]),ns)))
            
            % PHASE BY POSITION
            p(4, ns).select();
            if size(data.ThPrecess{i,ns},1)>1
                scatter([data.ThPrecess{i,ns}(:,1);...
                    data.ThPrecess{i,ns}(:,1)],...
                    [data.ThPrecess{i,ns}(:,2);...
                    data.ThPrecess{i,ns}(:,2)+360],15,'Filled','k');
            end
            axis tight;box off;axis off;axis square
            title(sprintf('Corr: %4.2f p: %4.2f',...
                data.measures(i,contains(data.varnames,["PhcircLinCorr","Phpval"]),ns)))
            % plot slope
            hold on
            plot_corr_line(data.ThPrecess{i,ns}(:,1),[],...
                data.measures(i,contains(data.varnames,["slopeCpU","phaseOffset"]),ns))
            
%             p(5, ns).select();
%             ax = gca;
%             ts=data_video_spk(data_video_spk(:,6)==0,1);
%             x=data_video_spk(data_video_spk(:,6)==0,2);
%             spkts=data_video_spk(data_video_spk(:,6)==1,1);
%             [~,I]=unique(ts);
%             ts=ts(I);
%             x=x(I);
%             Phasemap(ts,x,spkts,data.lfp,length(ratemap),trodeID,data.samplerate)
%             colormap(ax,jet(255))
            
            % PLOT AUTOCORRELATION
            p(5, ns).select();
            y=data.thetaautocorr{i,ns};
            plot(y,'LineWidth',2, 'color','k');
            axis tight
            hold on;box off; axis off
            title(sprintf('Theta modulation: %4.2f',...
                data.measures(i,contains(data.varnames,'thetaindex'),ns)))
            
            % PLOT AVERAGE WAVEFORMS [make this compatible with less channels]
            if ns==1
                p(6, ns).select();
                ax = gca;
                waves=zeros(4,150);
                waves(1:size(data.avgwave{i},1),:)=data.avgwave{i};
                plot(1:150,waves(1,:),'LineWidth',2, 'color','k');hold on
                plot(151:300,waves(2,:),'LineWidth',2, 'color','k');
                plot(1:150,waves(3,:)+max(abs([waves(1,:),waves(2,:)])),'LineWidth',2, 'color','k');
                plot(151:300,waves(4,:)+max(abs([waves(1,:),waves(2,:)])),'LineWidth',2, 'color','k');
                axis tight
                hold on;box off; axis off
                title(sprintf('ShortISI: %4.2f :IsoDist %4.2f',...
                    data.measures(i,contains(data.varnames,["IsolationDistance","ShortISI"]),ns)))
            end
        end
    end
%     F=getframe(gcf);
%     [X(:,:,:,i),~]=frame2im(F);
end

% summary plot
if plotraster==1
    % plot raster
    fig=figure('Name',[data.rat,'  ',data.sessionID],'NumberTitle','off');
    fig.Color=[1 1 1];
    subplot(3,2,[1 2])
    X = cellfun(@transpose,data.Spikes,'un',0);
    [~,I] = sort(cellfun(@length,X));
    X = X(I);
    plotSpikeRaster(X,'PlotType','vertline');hold on
    xlabel('Time(sec)')
    ylabel('Cells')
    title('Raster')
    for i=1:size(data.events,2)
        plot([data.events(1,i),data.events(1,i)],ylim,'b','LineWidth',2)
        plot([data.events(2,i),data.events(2,i)],ylim,'r','LineWidth',2)
    end
    % plot histogram
    subplot(3,2,[3 4])
    spiketimes=vertcat(data.Spikes{:});
    % BIN SPIKE TIMES FROM ALL CELLS OVER 500MS
    edges=linspace(min(spiketimes),max(spiketimes),round(data.frames(end,1)/3));
%     plot(edges(1:end-1),histcounts(spiketimes,edges),'k');
    histogram(spiketimes,edges)
    box off;axis tight
    xlabel('Time(sec)')
    ylabel('Count')
    title('All Cells Binned over 3 sec')    
    
    subplot(3,2,5)
    % find rows with all zeros
    idx=sum(data.lfp.signal==0,2)~=size(data.lfp.signal,2);
    % power spectrum
    pspectrum(median(data.lfp.signal(idx,:)),1000,'spectrogram','FrequencyLimits',[0, 100],'TimeResolution', 10);
    
    subplot(3,2,6)
    % power by freq
    pspectrum(median(data.lfp.signal(idx,:)),1000,'power','FrequencyLimits',[0, 100]);
end
% if xcorr
%     fig=figure('Name',[data.rat,'  ',data.sessionID,'_xcorr'],'NumberTitle','off');
%     fig.Color=[1 1 1];
%     cross_corr_matrix(data.Spikes,data.spikesID)
%     disp('Figures Done')
% end
end

%% SUB FUNCTIONS %%
function plot_corr_line(x,maxx,data)
if isstruct(data)
    for f=1:length(data.fields)
        % restrict x to field
        width=rescale([1:maxx,data.fields{f}.start,data.fields{f}.stop],0,1);
        width=width(end-1:end);
        xtemp=x(x>width(1) & x<width(2));
        
        % find x dims
        dx=max(xtemp)-min(xtemp);
        xm=mean([min(xtemp);max(xtemp)]);
        
        % kempter (kempter method usually finds a better fit than FMA)
        ym=wrapTo2Pi((data.fields{f}.ThPrecess.slopeCpU*2*pi)*xm...
            +data.fields{f}.ThPrecess.phaseOffset)*180/pi;
        X = xm+[-dx dx]/2;
        Y = ym+[-dx dx]/2*((data.fields{f}.ThPrecess.slopeCpU*2*pi)*180/pi);
        plot(X,Y,'r','linewidth',1.5);
        plot(X,Y+360,'r','linewidth',1.5);
        ylim([0 (4*pi)*180/pi])
        
        % FMA
%         if isfield(data.fields{f}.ThPrecess.stats,'slope')
%             ym = wrapTo2Pi(data.fields{f}.ThPrecess.stats.slope*xm...
%                 +data.fields{f}.ThPrecess.stats.intercept)*180/pi;
%             X = xm+[-dx dx]/2;
%             Y = ym+[-dx dx]/2*(data.fields{f}.ThPrecess.stats.slope*180/pi);
%             plot(X,Y,'r','linewidth',2);
%             plot(X,Y+360,'r','linewidth',2);
%             ylim([0 (4*pi)*180/pi])
%         end
    end
else
    % for open field 2d data using the pass index method
    % find x dims
    dx=max(x)-min(x);
    xm=mean([min(x);max(x)]);
    % kempter
    % y dims
    ym=wrapTo2Pi((data(1)*2*pi)*xm...
        +data(2))*180/pi;
    % find line
    X = xm+[-dx dx]/2;
    Y = ym+[-dx dx]/2*((data(1)*2*pi)*180/pi);
    plot(X,Y,'r','linewidth',2);
    plot(X,Y+360,'r','linewidth',2);
    ylim([0 (4*pi)*180/pi])
end
 end
 
 
%         end
        % plot filled and smoothed rate map
%         imAlpha=ones(size([SmoothRateMap;SmoothRateMap]));
%         imAlpha(isnan([SmoothRateMap;SmoothRateMap]))=0;
%         subplot(4,1,2); imagesc([SmoothRateMap;SmoothRateMap],'AlphaData',imAlpha);
%         hold on; shading interp; colormap jet(255); axis off; box off;axis image
%         title(['InfoContent: ',num2str(round(InformationContent,3)),' F2W: ',num2str(round(Field2Wall,3))]);
%         
%         subplot(4,1,3), area(SmoothRateMap(1,:),'LineWidth',2,'EdgeColor',[0,0,0]+0.4,'FaceColor',[0,0,0]+0.8);
%         box off; xlim([1 nBinsx]); hold on
%         % plot field location
%         if sum(field)>0
%             FL=find(field==1);
%             plot([FL(1);FL(1)],[PeakRate;0],'LineWidth', 2, 'color','r')
%             plot([FL(end);FL(end)],[PeakRate;0],'LineWidth', 2, 'color','r')
%         end
%         title(['PR: ',num2str(round(PeakRate,3)),' OFR: ',num2str(round(OverallFR,3))]);
%         set(fig,'Position',[842 42 838 924]);
%         
%         % PLOT SCATTERED PHASE PRECESSION
%         if ~isnan(ThPrecess.scatteredPH)
%             subplot(4,1,4)
%             % rescale x distance for plotting
%             %                         fieldidx=find(field_restrict);
%             if ~isempty(fieldidx)
%                 x=rescale(ThPrecess.scatteredPH(:,1),min(fieldidx),max(fieldidx));
%             else
%                 x=NaN;
%             end
%             plot(x,ThPrecess.scatteredPH(:,2),'k.');hold on;
%             plot(x,ThPrecess.scatteredPH(:,2)+360,'r.')
%             ylim([0 720]); xlim([1 length(SmoothRateMap)])
%             set(gca, 'YTick', [0;240;480;720],'Box','off');
%             title(['Precess R^2: ',num2str(round(ThPrecess.RSquared,3)),', Slope: ',...
%                 num2str(round(ThPrecess.slope,3))])
%             
%             % PLOT SMOOTHED RATEMAP
%             %                         subplot(5,1,5);
%             %                         h=pcolor(ThPrecess.smoothedPHmap);
%             %                         shading interp; colormap jet(255); hold on; box off; axis off; set(h, 'EdgeColor', 'none');
%             %                         title(['DOM: ',num2str(round(ThPrecess.DOM,3))])
%             
%             % PLOT AUTOCORR
%             %                         plot([linspace(-500,-1,100),0,linspace(1,500,100)], cor,'k'); hold on
%             %
%             %                         set(gca,'YTickLabel',[],'YTick',[],'XMinorTick','on','YMinorTick','off','LineWidth',1)
%             %                         line([0 0], ylim, 'linestyle', ':', 'color', [.7 .7 .7]);
%             % %                         set(gca, 'fontsize', 20, 'box', 'off');
%             %                         title(['Theta Index= ',num2str(thetaindex(j)),' Freq= ',num2str(peak)])
%         end
%         hold off
%         
%         % FOR CIRCULAR ARENA
%         %         elseif isequal(linear_track,'no')
%         % plot spike on path
%         frames=[rescale(data_video_spk(:,2),1,length(SmoothRateMap)),rescale(data_video_spk(:,3),1,length(SmoothRateMap))];
%         spkfram=frames(data_video_spk(:,6)==1,:);
%         fig=figure;subplot(3,2,1);
%         plot(frames(:,1), frames(:,2),'LineWidth',1,'color','k');
%         hold on; axis off
%         scatter(spkfram(:,1), spkfram(:,2), 35, 'filled', 'r');
%         box off; axis image
%         plot(bound(:,1),bound(:,2),'LineWidth', 3, 'color', [.3 .3 .3])
%         title([cell{end},' Cell: ',num2str(i),'nSpikes: ',...
%             num2str(nSpikes),' F2W: ',num2str(round(Field2Wall,3))]);
%         
%         imAlpha=ones(size(SmoothRateMap));
%         imAlpha(isnan(SmoothRateMap))=0;
%         subplot(3,2,2);imagesc(SmoothRateMap,'AlphaData',imAlpha);
%         axis xy; colormap jet(255); axis off; hold on; box off; axis image;
%         title(['IC: ',num2str(round(InformationContent,3)),' PR: ',num2str(round(PeakRate,3)),' OFR: ',num2str(round(OverallFR,3))]);
%         
%         %smooth spike data
%         [pdf,~]=circ_ksdensity(spks_VEL(:,4), 0:359,'msni');
%         % Plot Smoothed firing rate x HEAD DIRECTION
%         subplot(3,2,3); plot(pdf,'LineWidth',2,'color','k')
%         axis tight; hold on; xlim ([0 360]); box off
%         title(['PR: ',num2str(round(peak_Firing_Rate,3)),' D IC: ',num2str(round(Direct_infoContent,3))])
%         
%         % Firing rate x HD polar plot for the nonsmoothed data above
%         subplot(3,2,4); polarplot = polar(BinsAngle([1:60 1]),BinsNbSpikes([1:60 1]),'b');
%         set(polarplot,'linewidth',3,'color','k'); axis off
%         title(['RLength: ',num2str(round(mean_vector_length,3)),' Pref Dir: ',num2str(round(preferred_Direction,3))]);
%         set(0,'Showhiddenhandles','on')
%         set(fig,'Position',[842 42 838 924]);
%         %                     fig1 = figure(Session);
%         
%         if ~isnan(ThPrecess.scatteredPH)
%             subplot(3,2,5); plot(ThPrecess.scatteredPH(:,1),ThPrecess.scatteredPH(:,2),'k.');hold on;
%             plot(ThPrecess.scatteredPH(:,1),ThPrecess.scatteredPH(:,2)+360,'r.')
%             ylim([0 720]); xlim([0 1])
%             set(gca, 'YTick', [0;240;480;720],'Box','off');
%             title(['Precess',' R^2: ',num2str(round(ThPrecess.RSquared,3)),', Slope: ',...
%                 num2str(round(ThPrecess.slope,3)),', Corr: ',num2str(round(ThPrecess.lapCorrelation,3))])
%             hold off
%             
%             % PLOT SMOOTHED RATEMAP
%             %                         subplot(3,2,6); h=pcolor(ThPrecess.smoothedPHmap);
%             %                         shading interp; colormap jet(255); hold on; box off; axis off; set(h, 'EdgeColor', 'none');
%             %                         title(['DOM: ',num2str(round(ThPrecess.DOM,3))])
%         end
%         hold off
%         %         end
%         pause(.0000001)
% %     end
% % end
