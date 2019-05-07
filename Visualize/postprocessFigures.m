% postprocessFigures
% plots figures for every cell from standard bclarklab data structure with
% a call to of ( postprocessFigures.main(data) )
%
% You can also call any other of the listed methods shown below independently
%
% List of Methods
%   main
%   raster
%   plot_corr_line
%   unlinearpath
%   spikesonpath
%   spikesonpath_2d
%   ratemaps
%   ratemaps_2d
%   phase_by_pos
%   phase_by_pos_2d
%   phase_map
%   autocors
%   avg_waveforms
%   phase_colormap
%   plot_HD_tuning
%

% Ryan Harvey 2019

classdef postprocessFigures
    
    methods(Static)
        function main(data,cellid)

            % how many cells, how many sessions?
            ncells=length(data.Spikes);
            nsessions=size(data.events,2);
            if any(any(contains(data.mazetypes,'circ track') | contains(data.mazetypes,'LinearTrack')))
                nsessions=size(data.events,2)+1;
            end
            
            if sum(contains(data.mazetypes,'track'))>1
                nsessions=sum(contains(data.mazetypes,'track'))*2+sum(~contains(data.mazetypes,'track'));
            end
            % set up colormap
            %             theta=0:.01:2*pi;
            %             c=hsv(length(theta));
            
            %             if any(any(contains(data.mazetypes,'circ track') | contains(data.mazetypes,'LinearTrack')))
            %                 unlinear=data.linear_track.nonlinearFrames;
            %             end
            
            if exist('cellid','var')
                cells_to_find=strcat(cellid{1},num2str(cellid{2}));

                cell_list=strcat(data.spikesID.TetrodeNum,num2str(data.spikesID.CellNum));

                cells=find(ismember(cell_list,cells_to_find))';
                
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
                pause(.0001)
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
                    
                    if contains(data.mazetypes{round(ns/2)},'track','IgnoreCase',true)
                        
                        if ns==1 || ns==3% FOR LEFT RUNNING DIRECTION
                            postprocessFigures.unlinearpath(ax,data.linear_track{round(ns/2)}.nonlinearFrames,data.linear_track{round(ns/2)}.right{1,i}.laps,...
                                data.linear_track{round(ns/2)}.right{1,i}.dataspks,data.mazetypes,...
                                data.lfp.ts,data.lfp.theta_phase(trodeID,:))
                            
                            % PLOT SPIKES ON LINEARIZED PATH WITH EACH SPIKE COLOR
                            % CODED BY ITS THETA PHASE
                            p(2, ns).select();
                            ax=gca;
                            postprocessFigures.spikesonpath(ax,data.linear_track{round(ns/2)}.right{1,i}.dataspks,...
                                data.ratemap{i,ns},data.lfp.ts,data.lfp.theta_phase(trodeID,:))
                            
                            % PLOT EACH LAPS RATEMAP WITH THE OVERALL FIRING RATE
                            % SUPERIMPOSED
                            p(3, ns).select();
                            ax = gca;
                            postprocessFigures.ratemaps(ax,data.linear_track{round(ns/2)}.right{1,i}.maps,data.ratemap{i,ns},...
                                data.measures(i,:,ns),data.varnames)
                            
                            % PHASE BY POSITION
                            p(4, ns).select();
                            ax = gca;
                            postprocessFigures.phase_by_pos(ax,data.linear_track{round(ns/2)}.right{1,i}.dataspks,...
                                data.lfp.ts,data.lfp.theta_phase(trodeID,:),...
                                data.ratemap{i,ns},data.linear_track{round(ns/2)}.right{i})
                            
                            % PHASE MAP
                            p(5, ns).select();
                            ax = gca;
                            postprocessFigures.phase_map(ax,data.linear_track{round(ns/2)}.right{1,i}.dataspks(data.linear_track{round(ns/2)}.right{1,i}.dataspks(:,6)==0,1),...
                                data.linear_track{round(ns/2)}.right{1,i}.dataspks(data.linear_track{round(ns/2)}.right{1,i}.dataspks(:,6)==0,2),...
                                data.linear_track{round(ns/2)}.right{1,i}.dataspks(data.linear_track{round(ns/2)}.right{1,i}.dataspks(:,6)==1,1),...
                                data.lfp,data.ratemap{i,ns},trodeID,data.samplerate)
                            
                            % PLOT AUTOCORRELATION
                            p(6, ns).select();
                            postprocessFigures.autocors(data.thetaautocorr{i,ns},data.measures(i,:,ns),data.varnames)
                            
                            % PLOT AVERAGE WAVEFORMS
                            p(7, ns).select();
                            postprocessFigures.avg_waveforms(data.avgwave{i},data.measures(i,:,ns),data.varnames)
                            
                        elseif ns==2 || ns==4% FOR RIGHT RUNNING DIRECTION
                            
                            postprocessFigures.unlinearpath(ax,data.linear_track{round(ns/2)}.nonlinearFrames,data.linear_track{round(ns/2)}.left{1,i}.laps,...
                                data.linear_track{round(ns/2)}.left{1,i}.dataspks,data.mazetypes,...
                                data.lfp.ts,data.lfp.theta_phase(trodeID,:))
                            
                            % PLOT SPIKES ON LINEARIZED PATH WITH EACH SPIKE COLOR
                            % CODED BY ITS THETA PHASE
                            p(2, ns).select();
                            ax = gca;
                            postprocessFigures.spikesonpath(ax,data.linear_track{round(ns/2)}.left{1,i}.dataspks,...
                                data.ratemap{i,ns},data.lfp.ts,data.lfp.theta_phase(trodeID,:))
                            
                            % PLOT EACH LAPS RATEMAP WITH THE OVERALL FIRING RATE
                            % SUPERIMPOSED
                            p(3, ns).select();
                            ax = gca;
                            postprocessFigures.ratemaps(ax,data.linear_track{round(ns/2)}.left{1,i}.maps,data.ratemap{i,ns},...
                                data.measures(i,:,ns),data.varnames)
                            
                            % PHASE BY POSITION
                            p(4, ns).select();
                            ax = gca;
                            postprocessFigures.phase_by_pos(ax,data.linear_track{round(ns/2)}.left{1,i}.dataspks,...
                                data.lfp.ts,data.lfp.theta_phase(trodeID,:),...
                                data.ratemap{i,ns},data.linear_track{round(ns/2)}.left{i})
                            
                            % PHASE MAP
                            p(5, ns).select();
                            ax = gca;
                            postprocessFigures.phase_map(ax,data.linear_track{round(ns/2)}.left{1,i}.dataspks(data.linear_track{round(ns/2)}.left{1,i}.dataspks(:,6)==0,1),...
                                data.linear_track{round(ns/2)}.left{1,i}.dataspks(data.linear_track{round(ns/2)}.left{1,i}.dataspks(:,6)==0,2),...
                                data.linear_track{round(ns/2)}.left{1,i}.dataspks(data.linear_track{round(ns/2)}.left{1,i}.dataspks(:,6)==1,1),...
                                data.lfp,data.ratemap{i,ns},trodeID,data.samplerate)
                            
                            % PLOT AUTOCORRELATION
                            p(6, ns).select();
                            postprocessFigures.autocors(data.thetaautocorr{i,ns},data.measures(i,:,ns),data.varnames)
                            
                            % PHASE COLORMAP
                            p(7, ns).select();
                            ax=gca;
                            postprocessFigures.phase_colormap(ax)
                        end
                        
                    elseif (contains(data.mazetypes{round(ns/2)},'box','IgnoreCase',true) ||...
                            contains(data.mazetypes{round(ns/2)},'Cylinder','IgnoreCase',true)) &&...
                            nsessions>size(data.events,2) &&...
                            ns~=1 % FOR OPEN FIELD AND BOX
                        % PLOT TUNNING CURVE
                        p(1, ns).select();
                        [data_video_spk,~]=createframes_w_spikebinary(data,ns-1,i);
%                         postprocessFigures.plot_HD_tuning(data_video_spk,data.samplerate)
                        postprocessFigures.plot_HD_tuning(data,ns-1,i)

                        
                        % PLOT SPIKES ON PATH WITH EACH SPIKE COLOR
                        % CODED BY ITS THETA PHASE
                        p(2, ns).select();
                        ax=gca;
                        postprocessFigures.spikesonpath_2d(ax,data_video_spk,data.lfp.ts,data.lfp.theta_phase(trodeID,:))
                        
                        % PLOT RATEMAP
                        p(3, ns).select();
                        ax = gca;
                        postprocessFigures.ratemaps_2d(ax,data.ratemap{i,ns},data.measures(i,:,ns),data.varnames)
                        
                        % PHASE BY POSITION
                        p(4, ns).select();
                        postprocessFigures.phase_by_pos_2d(data.ThPrecess{i,ns},data.ratemap{i,ns},data.measures(i,:,ns),data.varnames)
                        
                        % PLOT AUTOCORRELATION
                        p(5, ns).select();
                        postprocessFigures.autocors(data.thetaautocorr{i,ns},data.measures(i,:,ns),data.varnames)
                        
                    elseif contains(data.mazetypes{round(ns/2)},'Cylinder') || contains(data.mazetypes{round(ns/2)},'box','IgnoreCase',true)% IF FIRST MAZE IS NOT A TRACK
                        
                        % PLOT TUNNING CURVE
                        p(1, ns).select();
                        [data_video_spk,~]=createframes_w_spikebinary(data,ns,i);
%                         postprocessFigures.plot_HD_tuning(data_video_spk,data.samplerate)
                        postprocessFigures.plot_HD_tuning(data,ns,i)
                        
                        % PLOT SPIKES ON PATH WITH EACH SPIKE COLOR
                        % CODED BY ITS THETA PHASE
                        p(2, ns).select();
                        ax=gca;
                        postprocessFigures.spikesonpath_2d(ax,data_video_spk,data.lfp.ts,data.lfp.theta_phase(trodeID,:))
                        
                        % PLOT RATEMAP
                        p(3, ns).select();
                        ax=gca;
                        postprocessFigures.ratemaps_2d(ax,data.ratemap{i,ns},data.measures(i,:,ns),data.varnames)
                        
                        % PHASE BY POSITION
                        p(4, ns).select();
                        postprocessFigures.phase_by_pos_2d(data.ThPrecess{i,ns},data.ratemap{i,ns},data.measures(i,:,ns),data.varnames)
                        
                        % PLOT AUTOCORRELATION
                        p(5, ns).select();
                        postprocessFigures.autocors(data.thetaautocorr{i,ns},data.measures(i,:,ns),data.varnames)
                        
                        % PLOT AVERAGE WAVEFORMS
                        if ns==1
                            p(6, ns).select();
                            postprocessFigures.avg_waveforms(data.avgwave{i},data.measures(i,:,ns),data.varnames)
                        end
                    end
                end
            end
            
            % summary plot
            if plotraster==1
                postprocessFigures.raster(data)
            end
            % if xcorr
            %     fig=figure('Name',[data.rat,'  ',data.sessionID,'_xcorr'],'NumberTitle','off');
            %     fig.Color=[1 1 1];
            %     cross_corr_matrix(data.Spikes,data.spikesID)
            %     disp('Figures Done')
            % end
        end
        
        function raster(data)
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
            
            colormap(viridis(255))
        end
        
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
        
        function unlinearpath(ax,unlinear,laps,spkframes,mazetypes,lfp_ts,theta_phase)
            % set up colormap
            theta=0:.01:2*pi;
            color=hsv(length(theta));
            
            % unlinearpath: plot raw pos
            for lap=1:length(laps)
                unlinear_{lap,1}=unlinear(ismember(unlinear(:,1),laps{lap}(:,1)),:);
            end
            
            spkts=spkframes(spkframes(:,6)==1,1);
            unlinear_together=vertcat(unlinear_{:});
            if contains(mazetypes{1},'LinearTrack')
                [unlinear_together(:,2),unlinear_together(:,3)]=...
                    ParameterizeToLinearTrack2(unlinear_together(:,2),unlinear_together(:,3));
            end
            plot(unlinear_together(:,2),unlinear_together(:,3),'.k');hold on
            [ts,idx]=unique(unlinear_together(:,1));
            scatter(interp1(ts,unlinear_together(idx,2),spkts),interp1(ts,unlinear_together(idx,3),spkts),10,...
                interp1(lfp_ts,theta_phase,spkts,'linear'),'filled');
            
            box off;axis image; axis off
            colormap(ax,color)
            title(['nSpikes: ',num2str(length(spkts))]);
        end
        
        function spikesonpath(ax,dataspks,ratemap,lfp_ts,theta_phase)
            % PLOT SPIKES ON LINEARIZED PATH WITH EACH SPIKE COLOR
            % CODED BY ITS THETA PHASE
            % set up colormap
            theta=0:.01:2*pi;
            color=hsv(length(theta));
            
            ts=dataspks(:,1);
            x=dataspks(:,2);
            x=rescale(x,1,length(ratemap));
            spkbinary=logical(dataspks(:,6));
            plot(x,ts,'.k');
            axis tight
            hold on;box off; axis off
            scatter(x(spkbinary),ts(spkbinary),20,...
                interp1(lfp_ts,...
                theta_phase,ts(spkbinary)','linear'),'filled');
            colormap(ax,color)
        end
        
        function spikesonpath_2d(ax,dataspks,lfp_ts,theta_phase)
            % PLOT SPIKES ON PATH WITH EACH SPIKE COLOR
            % CODED BY ITS THETA PHASE
            % set up colormap
            theta=0:.01:2*pi;
            color=hsv(length(theta));
            
            ts=dataspks(:,1);
            y=dataspks(:,3);
            x=dataspks(:,2);
            spkbinary=logical(dataspks(:,6));
            plot(x,y,'.k');
            axis image
            hold on;box off; axis off
            scatter(x(spkbinary),y(spkbinary),20,interp1(lfp_ts,theta_phase,ts(spkbinary)','linear'),'filled');
            colormap(ax,color)
            title(['nSpikes: ',num2str(sum(dataspks(:,6)==1))]);
        end
        
        function ratemaps(ax,lapmaps,ratemap,measures,varnames)
            % PLOT EACH LAPS RATEMAP WITH THE OVERALL FIRING RATE
            % SUPERIMPOSED
            laps=reshape([lapmaps{:}],[],length(lapmaps));
            imagesc(flipud(imrotate(laps,90)));hold on;set(gca,'Ydir','Normal')
            plot(rescale(ratemap,1,size(laps,2)),'LineWidth',2, 'color','w');
            set(gca,'XTickLabel',[]);
            axis tight
            hold on;box off;
            colormap(ax,viridis(255))
            title(sprintf('Info Content: %4.2f Peak Rate: %4.2f',...
                measures(contains(varnames,["InformationContent","PeakRate"]))))
        end
        
        function ratemaps_2d(ax,ratemap,measures,varnames)
            % ratemaps_2d
            imAlpha=ones(size(ratemap));
            imAlpha(isnan(ratemap))=0;
            imagesc(ratemap,'AlphaData',imAlpha);
            axis xy; axis off; hold on; box off; axis image;
            colormap(ax,viridis(255))
            title(sprintf('Info Content: %4.2f Peak Rate: %4.2f',...
                measures(contains(varnames,["InformationContent","PeakRate"]))))
        end
        
        function phase_by_pos(ax,dataspks,lfp_ts,theta_phase,ratemap,trackinfo)
            % PHASE BY POSITION
            x=rescale(dataspks(:,2),0,1);
            
            phase=interp1(lfp_ts,theta_phase,dataspks(logical(dataspks(:,6)),1)','linear');
            scatter([x(logical(dataspks(:,6)));x(logical(dataspks(:,6)))],[phase';phase'+2*pi]*180/pi,15,'Filled','k');
            box off;axis off
            xlim([min(x) max(x)])
            set(ax,'YTick',0:360:720,'YTickLabel',[{'0'},{'2\pi'},{'4\pi'}])
            
            if isfield(trackinfo,'fields')
                for f=1:length(trackinfo.fields)
                    rho(f)=trackinfo.fields{f}.ThPrecess.circLinCorr;
                    pval(f)=trackinfo.fields{f}.ThPrecess.pval;
                end
                a3=[rho;pval];
                title(sprintf(repmat(' rho:%4.2f p:%4.2f,',1,length(rho)),a3(:)'))
                clear rho pval
                hold on
                % plot correlation lines
                postprocessFigures.plot_corr_line(x(logical(dataspks(:,6))),length(ratemap),trackinfo)
            end
        end
        
        function phase_by_pos_2d(ThPrecess,ratemap,measures,varnames)
            % PHASE BY POSITION
            if size(ThPrecess,1)>1
                scatter([ThPrecess(:,1);ThPrecess(:,1)],[ThPrecess(:,2);ThPrecess(:,2)+360],10,'Filled','k');
            end
            axis tight;box off;axis off;axis square
            title(sprintf('Corr: %4.2f p: %4.2f',...
                measures(contains(varnames,["PhcircLinCorr","Phpval"]))))
            hold on
            postprocessFigures.plot_corr_line(ThPrecess(:,1),length(ratemap),...
                measures(contains(varnames,["slopeCpU","phaseOffset"])))
        end
        
        function phase_map(ax,ts,x,spkts,lfp,ratemap,trodeID,samplerate)
            % PHASE MAP
            [~,I]=unique(ts);
            ts=ts(I);
            x=x(I);
            bins=length(ratemap);
            
            xedge=linspace(min(x),max(x),bins+1);
            phaseedge=linspace(0,720,bins+1);
            
            xspk=interp1(ts,x,spkts);
            
            phase=interp1(lfp.ts,...
                lfp.theta_phase(trodeID,:),spkts','linear');
            spkmap=histcounts2([xspk;xspk],[phase';phase'+2*pi]*180/pi,xedge,phaseedge);
            
            phase=interp1(lfp.ts,...
                lfp.theta_phase(trodeID,:),ts','linear');
            occ=histcounts2([x;x],[phase';phase'+2*pi]*180/pi,xedge,phaseedge);
            occ=occ/samplerate;
            
            phasemap=spkmap./occ;
            
            phasemap(isnan(phasemap)) = 0;
            phasemap(isinf(phasemap)) = 0;
            
            % SMOOTH
            h=4;
            myfilter = fspecial('gaussian',[4 12]*h, h);
            phasemap = imfilter([phasemap,phasemap,phasemap],myfilter,'replicate');
            phasemap=phasemap(:,bins:(bins-1)*2);
            
            % plot
            pcolor(flipud(rot90(phasemap)));shading flat;box off;axis off;axis tight;
            colormap(ax,viridis(255))
        end
        
        function autocors(thetaautocorr,measures,varnames)
            % PLOT AUTOCORRELATION
            plot(thetaautocorr,'LineWidth',2, 'color','k');
            axis tight
            hold on;box off; axis off
            title(sprintf('Theta modulation: %4.2f',...
                measures(contains(varnames,'thetaindex'))))
        end
        
        function avg_waveforms(avgwave,measures,varnames)
            % PLOT AVERAGE WAVEFORMS
            waves=zeros(4,150);
            waves(1:size(avgwave,1),:)=avgwave;
            plot(1:150,waves(1,:),'LineWidth',2, 'color','k');hold on
            plot(151:300,waves(2,:),'LineWidth',2, 'color','k');
            plot(1:150,waves(3,:)+max(abs([waves(1,:),waves(2,:)])),'LineWidth',2, 'color','k');
            plot(151:300,waves(4,:)+max(abs([waves(1,:),waves(2,:)])),'LineWidth',2, 'color','k');
            axis tight
            hold on;box off; axis off
            title(sprintf('ShortISI: %4.2f :IsoDist %4.2f',...
                measures(contains(varnames,["IsolationDistance","ShortISI"]))))
        end
        
        function phase_colormap(ax)
            % set up colormap
            theta=0:.01:2*pi;
            color=hsv(length(theta));
            
            imagesc(theta);hold on
            colormap(ax,color)
            set(ax,'XTick',linspace(0,length(theta),2),'XTickLabel',...
                [{'0'} {'2\pi'}],'TickLength',[0,0],'YTick',[.5 1 1.5],'YTickLabel',[])
            x = 0:1:length(theta);
            y = rescale(gaussmf(x,[median([1,length(theta)/2]) length(theta)/2]),.5,1.5);
            plot(x,y,'w')
            axis square
        end
        
        function plot_HD_tuning(data,session,cell)
            % PLOT TUNNING CURVE
            [data_video_spk,~]=createframes_w_spikebinary(data,session,cell);

            [r,~,Ispk,~,~,tuning]=tuningcurve(data_video_spk(data_video_spk(:,6)==0,4),...
                data_video_spk(data_video_spk(:,6)==1,4),data.samplerate);
            %             angBins=0:6:360;
            %             bin_centers=movmedian(angBins,2);
            %             bin_centers(1)=[];
            %             polarplot(deg2rad(bin_centers),tuning,'k')
            %             ax=gca;
            %             set(ax,'RGrid','off','ThetaGrid','off','ThetaTick',[0 90],...
            %                 'ThetaTickLabels',[0 90],'RTick',max(tuning),'RTickLabel',max(tuning),'LineWidth',3);
            %             axis tight
            
            angBins=0:6:360;
            bin_centers=movmedian(angBins,2);
            bin_centers(1)=[];
            Polarplot = polar(deg2rad(bin_centers),tuning,'b');
            set(Polarplot,'linewidth',1,'color','k');
            axis off
            set(0,'Showhiddenhandles','on')
            extrastuff = setdiff(get(gca,'children'),Polarplot);
            delete(extrastuff)
            hold on
            horizontal=line([-max(tuning) max(tuning)],[0 0]); % for running max and min
            vertical=line([0 0],[-max(tuning) max(tuning)]);
            set(horizontal,'linewidth',2,'color',[.4 .4 .4]);
            set(vertical,'linewidth',2,'color',[.4 .4 .4]);
            axis image
            
            uistack(horizontal,'bottom')
            uistack(vertical,'bottom')
            
            title(sprintf('r: %4.2f DIC: %4.2f' ,[r,Ispk]))
            
            h=fill(get(Polarplot, 'XData'), get(Polarplot, 'YData'),...
                'k');
            set(h,'FaceAlpha',.5)

        end
    end
end