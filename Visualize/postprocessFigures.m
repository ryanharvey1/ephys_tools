% postprocessFigures
% plots figures for every cell from standard bclarklab data structure with
% a call to of ( postprocessFigures.main(data) )
%
% You can also call any other of the listed methods shown below independently
%
% List of Methods
%   main
%       input:
%               data: full postprocessed data structure
%               optional
%                   cellid: id of cell or cells you want to plot
%                   colorcode: spike color choice of 'HD', 'phase', or 'r'
%   raster(data)
%
%   unlinearpath(data,session,cell,direction,colorcode)
%
%   spikesonpath(data,session,cell,direction,colorcode)
%
%   spikesonpath_2d(data,session,cell,colorcode)
%
%   ratemaps(data,session,cell,direction)
%
%   ratemaps_2d(data,session,cell)
%
%   phase_by_pos(data,session,cell,direction)
%
%   phase_by_pos_2d(data,session,cell)
%
%   phase_map(data,session,cell,direction)
%
%   phasemap_2d(data,session,cell)
%
%   autocors(data,session,cell,direction)
%
%   avg_waveforms(data,session,cell)
%
%   phase_colormap(colorcode)
%
%   plot_HD_tuning(data,session,cell)


% Ryan Harvey 2019

classdef postprocessFigures
    
    methods(Static)
        function main(data,varargin)
            
            p = inputParser;
            p.addParameter('cellid',[]);
            p.addParameter('colorcode','HD');
            p.parse(varargin{:});
            cellid = p.Results.cellid;
            colorcode = p.Results.colorcode;
            
            % how many cells, how many sessions?
            ncells=length(data.Spikes);
            nsessions=size(data.events,2);
            
            sessions = 1:nsessions;
            for i = 2:nsessions
                if contains(data.mazetypes{i-1},'track','IgnoreCase',true)
                    sessions(i) = i+1;
                else
                    sessions(i) = sessions(i-1)+1;
                end
            end
        
            if any(contains(data.mazetypes,'track','IgnoreCase',true))
                nsessions=sum(contains(data.mazetypes,'track','IgnoreCase',true))*2+...
                    sum(~contains(data.mazetypes,'track','IgnoreCase',true));
            end
            
            if ~isempty(cellid)
                cells_to_find=strcat(cellid{1},num2str(cellid{2}));
                
                cell_list=strcat(data.spikesID.TetrodeNum,num2str(data.spikesID.CellNum));
                
                cells=find(ismember(strrep(cell_list,' ',''),cells_to_find))';
                
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
                % SET UP FIGURE
                pause(.0001)
                fig=figure('Name',[data.rat,'  ',data.sessionID,'  ',...
                    tetrode,' Cell: ',num2str(data.spikesID.CellNum(i))],...
                    'NumberTitle','off');
                fig.Color=[1 1 1];
                p = panel(fig);
                p.pack(7, nsessions);
                % set margins
                p.de.margin = 7;
                
                % and some properties
                p.fontsize = 12;
                p.fontweight='bold';
                
                true_num_sess = 1;
                
                for ns=sessions
                    
                    if contains(data.mazetypes{true_num_sess},'track','IgnoreCase',true)
                        direction = {'right','left'};
                        for dir_ = 1:2
                            % left running direction (labels are flipped)
                            p(1, ns+(dir_-1)).select();
                            postprocessFigures.unlinearpath(data,true_num_sess,i,direction{dir_},colorcode)
                            
                            % PLOT SPIKES ON LINEARIZED PATH WITH EACH SPIKE COLOR
                            % CODED BY ITS THETA PHASE
                            p(2, ns+(dir_-1)).select();
                            postprocessFigures.spikesonpath(data,true_num_sess,i,direction{dir_},colorcode)
                            
                            % PLOT EACH LAPS RATEMAP WITH THE OVERALL FIRING RATE
                            % SUPERIMPOSED
                            p(3, ns+(dir_-1)).select();
                            postprocessFigures.ratemaps(data,true_num_sess,i,direction{dir_})
                            
                            % PHASE BY POSITION
                            p(4, ns+(dir_-1)).select();
                            postprocessFigures.phase_by_pos(data,true_num_sess,i,direction{dir_})
                            
                            % PHASE MAP
                            p(5, ns+(dir_-1)).select();
                            postprocessFigures.phase_map(data,true_num_sess,i,direction{dir_})
                            
                            % PLOT AUTOCORRELATION
                            p(6, ns+(dir_-1)).select();
                            postprocessFigures.autocors(data,ns+(dir_-1),i)
                            
                            % PLOT AVERAGE WAVEFORMS
                            if ns == 1
                                p(7, ns).select();
                                postprocessFigures.avg_waveforms(data,true_num_sess,i)
                            end
                            
                            % PHASE COLORMAP
                            if dir_ == 2 && ns < 3
                                p(7, ns+(dir_-1)).select();
                                postprocessFigures.phase_colormap(colorcode)
                            end
                        end 
                        
                    else % else box or cylinder
                        
                        % PLOT TUNNING CURVE
                        p(1, ns).select();
                        postprocessFigures.plot_HD_tuning(data,true_num_sess,i)
                        
                        % PLOT SPIKES ON PATH WITH EACH SPIKE COLOR
                        % CODED BY ITS THETA PHASE
                        p(2, ns).select();
                        postprocessFigures.spikesonpath_2d(data,true_num_sess,i,colorcode)
                        
                        % PLOT RATEMAP
                        p(3, ns).select();
                        postprocessFigures.ratemaps_2d(data,ns,i)
                        
                        % PHASE BY POSITION
                        p(4, ns).select();
                        postprocessFigures.phase_by_pos_2d(data,ns,i)
                        
                        % plot phase map
                        p(5, ns).select();
                        postprocessFigures.phasemap_2d(data,true_num_sess,i)
                        
                        % PLOT AUTOCORRELATION
                        p(6, ns).select();
                        postprocessFigures.autocors(data,ns,i)
                        
                        % PLOT AVERAGE WAVEFORMS
                        if ns==1
                            p(7, ns).select();
                            postprocessFigures.avg_waveforms(data,true_num_sess,i)
                        end
                    end
                    true_num_sess = true_num_sess + 1;
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
            clear x y
            for i=1:length(X)
                [x{i},y{i}]=plotSpikeRaster(X(i),'PlotType','vertline');
            end
            for i=1:length(x)
                plot(x{i},y{i}+i,'Color',[rand(1),rand(1),rand(1)]);hold on
            end
            axis tight
            grid on
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
            % BIN SPIKE TIMES FROM ALL CELLS
            edges=min(spiketimes):3:max(spiketimes);
            yyaxis left
            h=histogram(spiketimes,edges);
            h.FaceColor = [.7 .7 .7];
            h.EdgeColor = [.7 .7 .7];
            box off;axis tight
            xlabel('Time(sec)')
            ylabel('Spike Count')
            title('All cells 3 sec bins')
            hold on
            yyaxis right
            plot(edges,smoothdata(interp1(data.frames(:,1),...
                data.frames(:,end),edges),'movmedian',40));
            ylabel('cm/sec')
            
            subplot(3,2,5)
            % find rows with all zeros
            [~,b,c]=unique(data.spikesID.TetrodeNum);
            idx=str2double(regexp(data.spikesID.TetrodeNum{b(mode(c))},'\d*','Match'));
            % power spectrum
            pspectrum(data.lfp.signal(idx,:),1000,'spectrogram','FrequencyLimits',...
                [0, 100],'TimeResolution', 10);
            
            subplot(3,2,6)
            % power by freq
            pspectrum((data.lfp.signal(idx,:)),1000,'power','FrequencyLimits',[0, 100]);
            
            colormap(viridis(255))
            darkBackground(gcf,[0.1 0.1 0.1],[0.7 0.7 0.7])
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
        
        function unlinearpath(data,session,cell,direction,colorcode)
            % unlinearpath: plot raw pos
            
            unlinear = data.linear_track{session}.nonlinearFrames;
            laps = data.linear_track{session}.(direction){1,cell}.laps;
            spkframes = data.linear_track{session}.(direction){1,cell}.dataspks;
            mazetypes = data.mazetypes;
            
            tetrode=strsplit(data.spikesID.paths{cell},filesep);
            tetrode=tetrode{end};
            trodeID=str2double(extractBetween(tetrode,'TT','.'));
            
            if ~iscell(laps) || isempty(laps)
                return
            end
            % set up colormap
            theta=0:.01:2*pi;
            color=hsv(length(theta));
            
            for lap=1:length(laps)
                unlinear_{lap,1}=unlinear(ismember(unlinear(:,1),laps{lap}(:,1)),:);
            end
            
            spkts=spkframes(spkframes(:,6)==1,1);
            unlinear_together=vertcat(unlinear_{:});
            if contains(mazetypes{1},'LinearTrack')
                [~,unlinear_together(:,2:3),~] = pca([unlinear_together(:,2),...
                    unlinear_together(:,3)]);
            end
            plot(unlinear_together(:,2),unlinear_together(:,3),'.k');hold on
            [ts,idx]=unique(unlinear_together(:,1));
            if strcmp(colorcode,'HD')
                scatter(interp1(ts,unlinear_together(idx,2),spkts),...
                    interp1(ts,unlinear_together(idx,3),spkts),10,...
                    interp1(rad2deg(theta)',color,...
                    interp1(ts,unlinear_together(idx,4),spkts)),'filled');
            elseif strcmp(colorcode,'phase')
                scatter(interp1(ts,unlinear_together(idx,2),spkts),...
                    interp1(ts,unlinear_together(idx,3),spkts),10,...
                    interp1(theta',color,...
                    interp1(data.lfp.ts,data.lfp.theta_phase(trodeID,:),spkts)),'filled');
            elseif strcmp(colorcode,'r')
                scatter(interp1(ts,unlinear_together(idx,2),spkts),...
                    interp1(ts,unlinear_together(idx,3),spkts),10,...
                    'r','filled');
            end
            box off;axis image; axis off
            title(['nSpikes: ',num2str(length(spkts))]);
        end
        
        function spikesonpath(data,session,cell,direction,colorcode)
            % PLOT SPIKES ON LINEARIZED PATH WITH EACH SPIKE COLOR
            % CODED BY ITS THETA PHASE
            
            dataspks = data.linear_track{session}.(direction){1,cell}.dataspks;
            
            ratemap = data.ratemap{cell,session};
            
            tetrode=strsplit(data.spikesID.paths{cell},filesep);
            tetrode=tetrode{end};
            trodeID=str2double(extractBetween(tetrode,'TT','.'));
            
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
            if strcmp(colorcode,'HD')
                scatter(x(spkbinary),ts(spkbinary),20,...
                    interp1(rad2deg(theta)',color,dataspks(spkbinary,4)),'filled');
            elseif strcmp(colorcode,'phase')
                scatter(x(spkbinary),ts(spkbinary),20,...
                    interp1(theta',...
                    color,interp1(data.lfp.ts,data.lfp.theta_phase(trodeID,:),...
                    ts(spkbinary)')),'filled');
            elseif strcmp(colorcode,'r')
                scatter(x(spkbinary),ts(spkbinary),20,'r','filled');
            end
        end
        
        function spikesonpath_2d(data,session,cell,colorcode)
            % PLOT SPIKES ON PATH WITH EACH SPIKE COLOR
            % CODED BY ITS THETA PHASE
            tetrode=strsplit(data.spikesID.paths{cell},filesep);
            tetrode=tetrode{end};
            trodeID=str2double(extractBetween(tetrode,'TT','.'));
            
            ax=gca;
            [dataspks,~] = createframes_w_spikebinary(data,session,cell);
            
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
            if strcmp(colorcode,'HD')
                scatter(x(spkbinary),y(spkbinary),20,...
                    interp1(rad2deg(theta)',color,dataspks(spkbinary,4)),'filled');
            elseif strcmp(colorcode,'phase')
                scatter(x(spkbinary),y(spkbinary),20,...
                    interp1(theta',...
                    color,interp1(data.lfp.ts,data.lfp.theta_phase(trodeID,:),...
                    ts(spkbinary)')),'filled');
            elseif strcmp(colorcode,'r')
                scatter(x(spkbinary),y(spkbinary),20,'r','filled');
            end
            title(['nSpikes: ',num2str(sum(dataspks(:,6)==1))]);
        end
        
        function ratemaps(data,session,cell,direction)
            % PLOT EACH LAPS RATEMAP WITH THE OVERALL FIRING RATE
            % SUPERIMPOSED
            
            ax = gca;
            lapmaps = data.linear_track{session}.(direction){1,cell}.maps;
            if ~iscell(lapmaps) || isempty(lapmaps)
                return
            end
            
            dataspk = data.linear_track{session}.(direction){1,cell}.dataspks;
            [SmoothRateMap,~,~,occ,~] = bindata(dataspk(dataspk(:,6)==0,:),...
                data.samplerate,dataspk(dataspk(:,6)==1,:),1,data.maze_size_cm(session));
             spatial_information = place_cell_analysis.SpatialInformation(...
                 'ratemap',SmoothRateMap,...
                 'occupancy',occ,'n_spikes',sum(dataspk(:,6)==1));
            laps=reshape([lapmaps{:}],[],length(lapmaps));
            imagesc(flipud(imrotate(laps,90)));hold on;set(gca,'Ydir','Normal')
            plot(rescale(nanmean(flipud(imrotate(laps,90))),1,size(laps,2)),...
                'LineWidth',2, 'color','w');
            set(gca,'XTickLabel',[]);
            axis tight
            hold on;box off;
            colormap(ax,viridis(255))
            title(sprintf('IC: %4.2f  %4.2f hz',spatial_information,max(SmoothRateMap)))
        end
        
        function ratemaps_2d(data,session,cell)
            % ratemaps_2d
            ax=gca;
            ratemap = data.ratemap{cell,session};
            measures = data.measures(cell,:,session);
            imAlpha=ones(size(ratemap));
            imAlpha(isnan(ratemap))=0;
            imagesc(ratemap,'AlphaData',imAlpha);
            axis xy; axis off; hold on; box off; axis image;
            colormap(ax,viridis(255))
            if ~isempty(data.varnames)
                title(sprintf('IC: %4.2f  %4.2f hz',...
                    measures(contains(data.varnames,["InformationContent","PeakRate"]))))
            end
        end
        
        function phase_by_pos(data,session,cell,direction)
            
            tetrode=strsplit(data.spikesID.paths{cell},filesep);
            tetrode=tetrode{end};
            trodeID=str2double(extractBetween(tetrode,'TT','.'));
            
            ax = gca;
            dataspks = data.linear_track{session}.(direction){1,cell}.dataspks;
            trackinfo = data.linear_track{session}.(direction){cell};
            
            if isempty(dataspks)
                return
            end
            
            [SmoothRateMap,~,~,~,~] = bindata(dataspks(dataspks(:,6)==0,:),...
                data.samplerate,dataspks(dataspks(:,6)==1,:),1,data.maze_size_cm(session));
            
            % PHASE BY POSITION
            x=rescale(dataspks(:,2),0,1);
            
            phase=interp1(data.lfp.ts,data.lfp.theta_phase(trodeID,:),...
                dataspks(logical(dataspks(:,6)),1)','linear');
            scatter([x(logical(dataspks(:,6)));x(logical(dataspks(:,6)))],...
                [phase';phase'+2*pi]*180/pi,15,'Filled','k');
            
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
                postprocessFigures.plot_corr_line(x(logical(dataspks(:,6))),...
                    length(SmoothRateMap),trackinfo)
            end
        end
        
        function phase_by_pos_2d(data,session,cell)
            % PHASE BY POSITION
            ThPrecess = data.ThPrecess{cell,session};
            
            if size(ThPrecess,1)>1
                scatter([ThPrecess(:,1);ThPrecess(:,1)],...
                    [ThPrecess(:,2);ThPrecess(:,2)+360],10,'Filled','k');
            end
            axis tight;box off;axis off;axis square
            title(sprintf('Corr: %4.2f p: %4.2f',...
                data.measures(cell,contains(data.varnames,["PhcircLinCorr","Phpval"]),session)))
            hold on
            postprocessFigures.plot_corr_line(ThPrecess(:,1),...
                length(data.ratemap{cell,session}),...
                data.measures(cell,contains(data.varnames,["slopeCpU","phaseOffset"]),session))
        end
        
        function phase_map(data,session,cell,direction)
            ax = gca;
            dataspks = data.linear_track{session}.(direction){1,cell}.dataspks;
            
            [SmoothRateMap,~,~,~,~] = bindata(dataspks(dataspks(:,6)==0,:),...
                data.samplerate,dataspks(dataspks(:,6)==1,:),1,data.maze_size_cm(session));
            
            ts = dataspks(dataspks(:,6) == 0,1);
            x = dataspks(dataspks(:,6) == 0,2);
            spkts = dataspks(dataspks(:,6) == 1,1);
            
            tetrode=strsplit(data.spikesID.paths{cell},filesep);
            tetrode=tetrode{end};
            trodeID=str2double(extractBetween(tetrode,'TT','.'));
            
            if isempty(x)
                return
            end
            % PHASE MAP
            [~,I]=unique(ts);
            ts=ts(I);
            x=x(I);
            bins=length(SmoothRateMap);
            
            xedge=linspace(min(x),max(x),bins+1);
            phaseedge=linspace(0,720,bins+1);
            
            xspk=interp1(ts,x,spkts);
            
            phase=interp1(data.lfp.ts,...
                data.lfp.theta_phase(trodeID,:),spkts','linear');
            spkmap=histcounts2([xspk;xspk],[phase';phase'+2*pi]*180/pi,xedge,phaseedge);
            
            phase=interp1(data.lfp.ts,...
                data.lfp.theta_phase(trodeID,:),ts','linear');
            occ=histcounts2([x;x],[phase';phase'+2*pi]*180/pi,xedge,phaseedge);
            occ=occ/data.samplerate;
            
            phasemap=spkmap./occ;
            
            phasemap(isnan(phasemap)) = 0;
            phasemap(isinf(phasemap)) = 0;
            
            % SMOOTH
            h=4;
            myfilter = fspecial('gaussian',[4 12]*h, h);
            phasemap = imfilter([phasemap,phasemap,phasemap],myfilter,'replicate');
            phasemap=phasemap(:,bins:(bins-1)*2);
            
            % plot
            imagesc(flipud(rot90(phasemap)));shading flat;box off;axis off;axis tight;
            colormap(ax,viridis(255))
        end
        
        function phasemap_2d(data,session,cell)
            
            tetrode=strsplit(data.spikesID.paths{cell},filesep);
            tetrode=tetrode{end};
            trodeID=str2double(extractBetween(tetrode,'TT','.'));
            
            [data_video_spk,data_video_nospk]=createframes_w_spikebinary(data,session,cell);
            
            [SmoothRateMap,~,~,~,~] = bindata(data_video_spk(data_video_spk(:,6)==0,:),...
                data.samplerate,data_video_spk(data_video_spk(:,6)==1,:),0,data.maze_size_cm(session));
            
            binside=mean([range(data_video_nospk(:,2))/length(SmoothRateMap),...
                range(data_video_nospk(:,3))/length(SmoothRateMap)]);
            
            results=pass_index(data_video_nospk(:,1),data_video_nospk(:,2:3),...
                data_video_spk(data_video_spk(:,6)==1,1),...
                [data.lfp.ts(data.lfp.ts>=data.events(1,session) &...
                data.lfp.ts<=data.events(2,session))]',...
                [data.lfp.signal(trodeID,data.lfp.ts>=data.events(1,session) &...
                data.lfp.ts<=data.events(2,session))]',...
                'plots',0,'method','place','binside',round(binside),...
                'sample_along','arc_length');
            
            bins=length(SmoothRateMap);
            xedge=linspace(-1,1,bins+1);
            phaseedge=linspace(0,720,bins*4);
            
            phase_spk = interp1(data.lfp.ts,data.lfp.theta_phase(trodeID,:),...
                data_video_spk(data_video_spk(:,6)==1,1));
            
            x_spk = interp1(results.ts,results.pass_index,...
                data_video_spk(data_video_spk(:,6)==1,1));
            
            spkmap=histcounts2([x_spk;x_spk],[phase_spk;phase_spk+2*pi]*180/pi,...
                xedge,phaseedge);
            
            phase_occ = interp1(data.lfp.ts,data.lfp.theta_phase(trodeID,:),...
                data_video_spk(data_video_spk(:,6)==0,1));
            
            x_occ = interp1(results.ts,results.pass_index,...
                data_video_spk(data_video_spk(:,6)==0,1));
            
            occ=histcounts2([x_occ;x_occ],[phase_occ;phase_occ+2*pi]*180/pi,...
                xedge,phaseedge);
            
            occ = occ/data.samplerate;
            
            phasemap=spkmap./occ;
            
            phasemap(isnan(phasemap)) = 0;
            phasemap(isinf(phasemap)) = 0;
            
            h=4;
            phase_size = size(phasemap,2);
            myfilter = fspecial('gaussian',[4 24]*h, h);
            phasemap = imfilter([phasemap,phasemap,phasemap],myfilter,'replicate');
            phasemap = phasemap(:,phase_size:phase_size*2);
            
            imagesc(flipud(rot90(phasemap)));shading flat;box off;axis off;axis tight;
            colormap(gca,viridis(255))
            axis square
        end
        
        
        function autocors(data,session,cell)
            % PLOT AUTOCORRELATION
            measures = data.measures(cell,:,session);
            plot(data.thetaautocorr{cell,session},'LineWidth',2, 'color','k');
            axis tight
            hold on;box off; axis off
            title(sprintf('Theta mod: %4.2f',...
                measures(contains(data.varnames,'thetaindex'))))
        end
        
        function avg_waveforms(data,session,cell)
            % PLOT AVERAGE WAVEFORMS
            
            avgwave = data.avgwave{cell};
            measures = data.measures(cell,:,session);
            
            sam=size(avgwave,2);
            waves=zeros(4,sam);
            waves(1:size(avgwave,1),:)=avgwave;
            plot(1:sam,waves(1,:),'LineWidth',2, 'color','k');hold on
            plot(sam+1:sam*2,waves(2,:),'LineWidth',2, 'color','k');
            plot(1:sam,waves(3,:)+max(abs([waves(1,:),waves(2,:)])),...
                'LineWidth',2, 'color','k');
            plot(sam+1:sam*2,waves(4,:)+max(abs([waves(1,:),waves(2,:)])),...
                'LineWidth',2, 'color','k');
            axis tight
            hold on;box off; axis off
            title(sprintf('ShortISI: %4.2f :IsoDist %4.2f',...
                measures(contains(data.varnames,["IsolationDistance","ShortISI"]))))
        end
        
        function phase_colormap(colorcode)
            % input:
            %   colorcode: just the title of the plot
            %
            % Set parameters (these could be arguments to a function)
            rInner = 250;     % inner radius of the colour ring
            rOuter = 500;    % outer radius of the colour ring
            ncols = 255;      % number of colour segments
            % Get polar coordinates of each point in the domain
            [x, y] = meshgrid(-rOuter:rOuter);
            [theta, rho] = cart2pol(x, y);
            % Set up colour wheel in hsv space
            hue = (theta + pi) / (2 * pi);     % hue into range (0, 1]
            hue = ceil(hue * ncols) / ncols;   % quantise hue
            saturation = ones(size(hue));      % full saturation
            brightness = double(rho >= rInner & rho <= rOuter);  % black outside ring
            % Convert to rgb space for display
            rgb = hsv2rgb(cat(3, hue, saturation, brightness));
            for i=1:3
                temp=rgb(:,:,i);
                temp(~brightness)=1;
                rgb(:,:,i)=temp;
            end
            imshow(flip(rgb,2));
            axis image
            if exist('colorcode','var')
                title(colorcode)
            end
        end
        
        function plot_HD_tuning(data,session,cell)
            % PLOT TUNNING CURVE
            [data_video_spk,~]=createframes_w_spikebinary(data,session,cell);
            
            [r,~,Ispk,~,~,tuning]=tuningcurve(data_video_spk(data_video_spk(:,6)==0,4),...
                data_video_spk(data_video_spk(:,6)==1,4),data.samplerate);
            
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