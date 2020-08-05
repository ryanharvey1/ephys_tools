classdef ripple_figs
    methods(Static)
        
        function lfp_viewer(ripple_info,data)
            
            if ~exist('data','var')
                disp('loading lfp ...')
                data = load(ripple_info.ripples.detectorinfo.ProcessedDatafile,'lfp');
            end
            figure
            trip = colormap(cool(size(data.lfp.signal,1)));
            for i = 1:size(data.lfp.signal,1)
                plot(data.lfp.ts,zscore(data.lfp.signal(i,:)) + i*2,'Color',trip(i,:))
                hold on
            end
            EMGFromLFP = load(ripple_info.ripples.detectorinfo.detectionparms.EMGfilename);
            plot(EMGFromLFP.timestamps,zscore(EMGFromLFP.data),'k')

            for r = 1:size(ripple_info.ripples.filtered_ripple,2)
                
                plot([ripple_info.ripples.peaks(r),ripple_info.ripples.peaks(r)],...
                    [min(ylim) max(ylim)],'r')
                hold on
                xbars = ripple_info.ripples.timestamps(r,:);
                patch([xbars(1) xbars(1), xbars(2) xbars(2)],...
                    [min(ylim) max(ylim) max(ylim) min(ylim)],...
                    [1 0 0],'FaceAlpha',.3,'EdgeColor','none')
            end
            
            ch_used = num2str(ripple_info.ripples.detectorinfo.detectionparms.channel_used);
            noise = num2str(ripple_info.ripples.detectorinfo.detectionparms.noise_channel);  
            title(['ch used: ',ch_used,' noise ch: ',noise,...
                ', Black: EMG',', Red: ripple'])
            
            darkBackground(gcf,[0.1 0.1 0.1],[0.7 0.7 0.7])
        end
        
        function ripple_grid(ripple_info)
            figure;
            p = ceil(sqrt(size(ripple_info.ripples.unfiltered_ripple,2)));
            colors = viridis(length(unique(ripple_info.ripples.peakNormedPower)));
            for r = 1:size(ripple_info.ripples.filtered_ripple,2)
                subplot(p,p,r)
                plot(ripple_info.maps.unfiltered_ripples(r,:),...
                    'color',interp1(unique(ripple_info.ripples.peakNormedPower),...
                    colors,ripple_info.ripples.peakNormedPower(r)))
                hold on
                box off
                axis off
            end
            sgtitle([num2str(r),' Ripples on Channel ',...
                num2str(ripple_info.ripples.detectorinfo.detectionparms.channel_used),...
                ', Color code: power'],'Color','w')
            darkBackground(gcf,[0.1 0.1 0.1],[0.7 0.7 0.7])
        end
        
        function ripple_per_channel(ripple_info,data)
            
            if ~exist('data','var')
                disp('loading lfp ...')
                data = load(ripple_info.ripples.detectorinfo.ProcessedDatafile,'lfp');
            end
            
            n_ripples = size(ripple_info.ripples.timestamps,1);
            colors = [rand(n_ripples,1),rand(n_ripples,1),rand(n_ripples,1)];
            
            duration = 1;
            figure
            for r = 1:n_ripples
                idx = data.lfp.ts >= ripple_info.ripples.timestamps(r,1) &...
                    data.lfp.ts <= ripple_info.ripples.timestamps(r,2);
                plot(duration:duration+sum(idx)-1,...
                    zscore(data.lfp.signal(:,idx),[],2)+[[1:size(data.lfp.signal,1)]*3]',...
                    'Color',colors(r,:))
                box off
                hold on
                duration = sum(idx) + duration;
            end
            title('Unfiltered ripples')
            darkBackground(gcf,[0.1 0.1 0.1],[0.7 0.7 0.7])
            
            % filtered plot
            disp('filtering lfp ...')
            for i =1:size(data.lfp.signal,1)
                signal_filtered(i,:) = BandpassFilter(data.lfp.signal(i,:), 1000,...
                    ripple_info.ripples.detectorinfo.detectionparms.passband);
            end
            duration = 1;
            figure
            for r = 1:n_ripples
                idx = data.lfp.ts >= ripple_info.ripples.timestamps(r,1) &...
                    data.lfp.ts <= ripple_info.ripples.timestamps(r,2);
                plot(duration:duration+sum(idx)-1,...
                    zscore(signal_filtered(:,idx),[],2)+[[1:size(signal_filtered,1)]*3]',...
                    'Color',colors(r,:))
                box off
                hold on
                duration = sum(idx) + duration;
            end
            title('Filtered ripples')
            darkBackground(gcf,[0.1 0.1 0.1],[0.7 0.7 0.7])
        end
        
        function ripple_location(ripple_info,data)
            
            if ~exist('data','var')
                disp('loading lfp ...')
                data = load(ripple_info.ripples.detectorinfo.ProcessedDatafile,...
                    'mazetypes','linear_track','events','frames');
            end
            figure;
            for i = 1:size(data.mazetypes,2)
                subplot(1,size(data.mazetypes,2),i)
                if contains(data.mazetypes{i},"track",'IgnoreCase',true)
                    frames = data.linear_track{1}.nonlinearFrames;
                    plot(frames(:,2),frames(:,3),'Color',[.7 .7 .7])
                    hold on
                    ripple_frames = interp1(frames(:,1),frames,ripple_info.ripples.peaks);
                    scatter(ripple_frames(:,2),ripple_frames(:,3),10,'r','Filled')
                    axis image
                else
                    idx = data.frames(:,1) > data.events(1,i) &...
                        data.frames(:,1) < data.events(2,i);
                    frames = data.frames(idx,:);
                    plot(frames(:,2),frames(:,3),'Color',[.7 .7 .7])
                    hold on
                    ripple_frames = interp1(frames(:,1),frames,ripple_info.ripples.peaks);
                    scatter(ripple_frames(:,2),ripple_frames(:,3),10,'r','Filled')
                    axis image
                end
                title(data.mazetypes{i})
            end
            sgtitle('Red dot = ripple','color','r')
            darkBackground(gcf,[0.1 0.1 0.1],[0.7 0.7 0.7])
        end
    end
end