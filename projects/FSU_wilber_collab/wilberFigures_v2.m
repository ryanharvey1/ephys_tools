function wilberFigures_v2(data)

if isfield(data,'notes')
    if contains(data.notes.trials,'rightleft')
        multidir(data)
    else
        singledir(data)
    end
else
    singledir(data)
end

    function singledir(data)
        for i=1:length(data.Spikes)
            fig=figure('Name',[data.rat,'  ',data.sessionID,'  ',' Cell: ',data.cellid{i}],'NumberTitle','off');
            fig.Color=[1 1 1];
            
            subplot(4,1,1)
            plot(data.frames_w_spk{i}(:,2),data.frames_w_spk{i}(:,1),'.k');hold on
            scatter(data.frames_w_spk{i}(data.frames_w_spk{i}(:,end)==1,2),...
                data.frames_w_spk{i}(data.frames_w_spk{i}(:,end)==1,1),200,'.r');
            box off;
            xlim([min(data.frames_w_spk{i}(:,2)) data.mazesize])
            xlabel('Track Length(cm)')
            ylabel('Time(sec)')
            title(sprintf('nSpikes: %4.0f',sum(data.frames_w_spk{i}(:,end))))
            xlim([0 120])
            
            subplot(4,1,2)
            ax = gca;
            imAlpha=ones(size(data.ratemaps{i}));
            imAlpha(isnan(data.ratemaps{i}))=0;
            imagesc(data.ratemaps{i},'AlphaData',imAlpha);
            axis xy; box off;
            xlabel('Track Length(bins)')
            ylabel('Trials')
            colormap(ax,jet(255))
            xlim([1 size(data.ratemaps{i},2)])
            title(sprintf('Peak Rate: %4.2fHz',max(data.ratemaps{i}(:))))
            
            
            subplot(4,1,3)
            if ~isempty(data.lfp)
                phase=interp1(data.lfp.ts,data.lfp.theta_phase,...
                    data.frames_w_spk{i}(data.frames_w_spk{i}(:,end)==1,1),'linear');
                scatter([data.frames_w_spk{i}(data.frames_w_spk{i}(:,end)==1,2);...
                    data.frames_w_spk{i}(data.frames_w_spk{i}(:,end)==1,2)],...
                    [phase',phase'+2*pi]*180/pi,15,'Filled','k');
                box off; axis tight
                xlim([min(data.frames_w_spk{i}(:,2)) data.mazesize])
                ylim([0 720])
                xlabel('Track Length(cm)')
                ylabel('Theta Phase')
            end
            subplot(4,1,4)
            plot(0:100,data.autocor(i,:),'k','LineWidth',2)
            set(gca,'XTick',0:10:100,'XTickLabel',[500:-100:100,0,100:100:500])
            axis xy; box off; axis tight
            xlabel('Lag(ms)')
            title(sprintf('Theta modulation: %4.2f',data.results.thetaindex(i)))
            
            
            set(findobj(gcf,'type','axes'),'FontSize',12,'FontWeight','Bold', 'LineWidth', 1.5);
        end
    end

    function multidir(data)
        for i=1:length(data.Spikes)
            fig=figure('Name',[data.rat,'  ',data.sessionID,'  ',' Cell: ',data.cellid{i}],'NumberTitle','off');
            fig.Color=[1 1 1];
            
            subplot(5,2,1.5)
            plot(data.frames_w_spk{i}(:,2),data.frames_w_spk{i}(:,1),'.k');hold on
            scatter(data.frames_w_spk{i}(data.frames_w_spk{i}(:,end)==1,2),...
                data.frames_w_spk{i}(data.frames_w_spk{i}(:,end)==1,1),200,'.r');
            box off;
            xlim([min(data.frames_w_spk{i}(:,2)) data.mazesize])
            xlabel('Track Length(cm)')
            ylabel('Time(sec)')
            title(sprintf('nSpikes: %4.0f',sum(data.frames_w_spk{i}(:,end))))
            xlim([0 120])
            
            subplot(5,2,3.5)
            ax = gca;
            imAlpha=ones(size(data.ratemaps{i}));
            imAlpha(isnan(data.ratemaps{i}))=0;
            imagesc(data.ratemaps{i},'AlphaData',imAlpha);
            axis xy; box off;
            xlabel('Track Length(bins)')
            ylabel('Trials')
            colormap(ax,jet(255))
            xlim([1 size(data.ratemaps{i},2)])
            title(sprintf('Peak Rate: %4.2fHz',max(data.ratemaps{i}(:))))
            
            subplot(5,2,5.5)
            plot(0:100,data.autocor(i,:),'k','LineWidth',2)
            set(gca,'XTick',0:10:100,'XTickLabel',[500:-100:100,0,100:100:500])
            axis xy; box off; axis tight
            xlabel('Lag(ms)')
            title(sprintf('Theta modulation: %4.2f',data.results.thetaindex(i)))
            
            subplot(5,2,7)
            plot(data.left{i}.frames_w_spk(:,2),data.left{i}.frames_w_spk(:,1),'.k');hold on
            scatter(data.left{i}.frames_w_spk(data.left{i}.frames_w_spk(:,end)==1,2),...
                data.left{i}.frames_w_spk(data.left{i}.frames_w_spk(:,end)==1,1),200,'.r');
            box off;
            xlim([min(data.left{i}.frames_w_spk(:,2)) data.mazesize])
            xlabel('Track Length(cm)')
            ylabel('Time(sec)')
            title(sprintf('left trials nSpikes: %4.0f',sum(data.left{i}.frames_w_spk(:,end))))
            xlim([0 120])
            
            subplot(5,2,8)
            plot(data.right{i}.frames_w_spk(:,2),data.right{i}.frames_w_spk(:,1),'.k');hold on
            scatter(data.right{i}.frames_w_spk(data.right{i}.frames_w_spk(:,end)==1,2),...
                data.right{i}.frames_w_spk(data.right{i}.frames_w_spk(:,end)==1,1),200,'.r');
            box off;
            xlim([min(data.right{i}.frames_w_spk(:,2)) data.mazesize])
            xlabel('Track Length(cm)')
            ylabel('Time(sec)')
            title(sprintf('right trials nSpikes: %4.0f',sum(data.right{i}.frames_w_spk(:,end))))
            xlim([0 120])
            
            subplot(5,2,9)
            ax = gca;
            imAlpha=ones(size(data.left{i}.ratemaps));
            imAlpha(isnan(data.left{i}.ratemaps))=0;
            imagesc(data.left{i}.ratemaps,'AlphaData',imAlpha);
            axis xy; box off;
            xlabel('Track Length(bins)')
            ylabel('Trials')
            colormap(ax,jet(255))
            xlim([1 size(data.left{i}.ratemaps,2)])
            title(sprintf('Peak Rate: %4.2fHz',max(data.left{i}.ratemaps(:))))
            
            subplot(5,2,10)
            ax = gca;
            imAlpha=ones(size(data.right{i}.ratemaps));
            imAlpha(isnan(data.right{i}.ratemaps))=0;
            imagesc(data.right{i}.ratemaps,'AlphaData',imAlpha);
            axis xy; box off;
            xlabel('Track Length(bins)')
            ylabel('Trials')
            colormap(ax,jet(255))
            xlim([1 size(data.right{i}.ratemaps,2)])
            title(sprintf('Peak Rate: %4.2fHz',max(data.right{i}.ratemaps(:))))
            
            set(findobj(gcf,'type','axes'),'FontSize',8,'FontWeight','Bold', 'LineWidth', 1.5);
        end
    end
end