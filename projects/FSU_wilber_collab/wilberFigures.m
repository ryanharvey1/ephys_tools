function wilberFigures(data)

% how many cells, how many sessions?
ncells=length(data.Spikes);
nsessions=size(data.events,2);
if any(any(contains(data.mazetype,'CircularTrack') | contains(data.mazetype,'LinearTrack')))
    nsessions=size(data.events,2)+1;
end
trodeID=1;
% set up colormap
theta=0:.01:2*pi;
c=hsv(length(theta));

if any(any(contains(data.mazetype,'CircularTrack') | contains(data.mazetype,'LinearTrack')))
    unlinear=data.linear_track.nonlinearFrames;
end

for i=1:ncells
    
    % SET UP FIGURE
    fig=figure('Name',[data.rat,'  ',data.sessionID,'  ',' Cell: ',data.cellid{i}],'NumberTitle','off');
    fig.Color=[1 1 1];
    p = panel(fig);
    p.pack(8, nsessions);
    % set margins
    p.de.margin = 7;
    
    % and some properties
    p.fontsize = 12;
    p.fontweight='bold';
    
    for ns=1:nsessions

        
        if ns==1 % FOR RIGHT RUNNING DIRECTION
            p(1, ns).select();
            ax=gca;
            % overall
            plot(data.frames(:,2),data.frames(:,3),'.k');hold on
            spkframes=data.frames_w_spk{i};
            spkts=spkframes(spkframes(:,end)==1,1);
            [ts,idx]=unique(data.frames(:,1));
            
            scatter(interp1(ts,data.frames(idx,2),spkts),interp1(ts,data.frames(idx,3),spkts),10,...
                interp1(data.lfp.ts,...
                data.lfp.theta_phase(trodeID,:),spkts,'linear'),'filled');
            box off;axis image; axis off
            colormap(ax,c)
            title(['Overall Trials, nSpikes: ',num2str(length(spkts))]);
            
            % overall ratemap
            p(2, ns).select();
            ax=gca;
            plot(data.overall_ratemap(i,:),'LineWidth',2, 'color','k');
            set(gca,'XTickLabel',[]);
            axis tight
            hold on;box off;
            title(sprintf('Info Content: %4.2f Peak Rate: %4.2f',...
                data.results.overall_InformationContent(i),max(data.overall_ratemap(i,:))))
            
            % raw path with directional spikes
            p(3, ns).select();
            ax=gca;
            spkframes=data.linear_track.right{1,i}.dataspks;
            spkts=spkframes(spkframes(:,end)==1,1);
            
            plot(data.frames(:,2),data.frames(:,3),'.k');hold on
            [ts,idx]=unique(data.frames(:,1));
            scatter(interp1(ts,data.frames(idx,2),spkts),interp1(ts,data.frames(idx,3),spkts),10,...
                interp1(data.lfp.ts,...
                data.lfp.theta_phase(trodeID,:),spkts,'linear'),'filled');
            box off;axis image; axis off
            colormap(ax,c)
            title(['Left Trials, nSpikes: ',num2str(length(spkts))]);
            
            % PLOT SPIKES ON LINEARIZED PATH WITH EACH SPIKE COLOR
            % CODED BY ITS THETA PHASE
            p(4, ns).select();
            ax=gca;
            ts=data.linear_track.right{1,i}.dataspks(:,1);
            x=data.linear_track.right{1,i}.dataspks(:,2);
            x=rescale(x,1,length(data.ratemap(ns,:,i)));
            spkbinary=logical(data.linear_track.right{1,i}.dataspks(:,end));
            plot(x,ts,'.k');
            axis tight
            hold on;box off; axis off
            scatter(x(spkbinary),ts(spkbinary),20,...
                interp1(data.lfp.ts,...
                data.lfp.theta_phase(trodeID,:),ts(spkbinary)','linear'),'filled');
            colormap(ax,c)
            
            % PLOT EACH LAPS RATEMAP WITH THE OVERALL FIRING RATE
            % SUPERIMPOSED
            p(5, ns).select();
            ax = gca;
            laps=reshape([data.linear_track.right{1,i}.maps{:}],...
                [],length(data.linear_track.right{1,i}.maps));
            imagesc(flipud(imrotate(laps,90)));hold on;set(gca,'Ydir','Normal')
            ratemap=data.ratemap(ns,:,i);
            plot(rescale(ratemap,1,size(laps,2)),'LineWidth',2, 'color','w');
            set(gca,'XTickLabel',[]);
            axis tight
            hold on;box off;
            colormap(ax,jet(255))
            title(sprintf('Info Content: %4.2f Peak Rate: %4.2f',...
                data.results.InformationContent(ns,i),data.results.PeakRate(ns,i)))
            
            
            % PHASE BY POSITION
            p(6, ns).select();
            ax = gca;
            x=rescale(x,0,1);
            phase=interp1(data.lfp.ts,...
                data.lfp.theta_phase(trodeID,:),ts(spkbinary)','linear');
            scatter([x(spkbinary);x(spkbinary)],[phase';phase'+2*pi]*180/pi,15,'Filled','k');
            axis tight;box off;axis off
            xlim([min(x) max(x)])
            
            for f=1:length(data.linear_track.right{1,i}.fields)
                rho(f)=data.linear_track.right{1,i}.fields{f}.ThPrecess.circLinCorr;
                pval(f)=data.linear_track.right{1,i}.fields{f}.ThPrecess.pval;
            end
            a3=[rho;pval];
            title(sprintf(repmat(' rho:%4.2f p:%4.2f,',1,length(rho)),a3(:)'))
            clear rho pval
            hold on
            % plot correlation lines
            %                 plot_corr_line(x(spkbinary),length(ratemap),data.linear_track.right{1,i})
            
            p(7, ns).select();
            ax = gca;
            ts=data.linear_track.right{1,i}.dataspks(data.linear_track.right{1,i}.dataspks(:,end)==0,1);
            x=data.linear_track.right{1,i}.dataspks(data.linear_track.right{1,i}.dataspks(:,end)==0,2);
            spkts=data.linear_track.right{1,i}.dataspks(data.linear_track.right{1,i}.dataspks(:,end)==1,1);
            [~,I]=unique(ts);
            ts=ts(I);
            x=x(I);
            Phasemap(ts,x,spkts,data.lfp,length(ratemap),trodeID,data.samplerate)
            colormap(ax,jet(255))
            
            % PLOT AUTOCORRELATION
            p(8, ns).select();
            ax = gca;
            y=data.autocor(ns,:,i);
            plot(y,'LineWidth',2, 'color','k');
            axis tight
            hold on;box off; axis off
            title(sprintf('Theta modulation: %4.2f',data.results.thetaindex(ns,i)))
            
        elseif ns==2 % FOR LEFT RUNNING DIRECTION
            p(1, ns).select();
            ax=gca;
            cb=imagesc(theta);hold on
            colormap(ax,c)
            set(ax,'XTick',linspace(0,length(theta),3),'XTickLabel',...
                [0 3.14 6.28],'TickLength',[0,0],'YTick',[.5 1 1.5],'YTickLabel',[])
            x = 0:1:length(theta);
            y = rescale(gaussmf(x,[median([1,length(theta)/2]) length(theta)/2]),.5,1.5);
            plot(x,y,'w')
            axis tight
            
            p(3, ns).select();
            ax=gca;
            % raw path
            spkframes=data.linear_track.left{1,i}.dataspks;
            spkts=spkframes(spkframes(:,end)==1,1);
            
            plot(data.frames(:,2),data.frames(:,3),'.k');hold on
            [ts,idx]=unique(data.frames(:,1));
            scatter(interp1(ts,data.frames(idx,2),spkts),interp1(ts,data.frames(idx,3),spkts),10,...
                interp1(data.lfp.ts,...
                data.lfp.theta_phase(trodeID,:),spkts,'linear'),'filled');
            box off;axis image; axis off
            colormap(ax,c)
            title(['right Trials, nSpikes: ',num2str(length(spkts))]);
            
            
            % PLOT SPIKES ON LINEARIZED PATH WITH EACH SPIKE COLOR
            % CODED BY ITS THETA PHASE
            p(4, ns).select();
            ax = gca;
            ts=data.linear_track.left{1,i}.dataspks(:,1);
            x=data.linear_track.left{1,i}.dataspks(:,2);
            x=rescale(x,1,length(data.ratemap(ns,:,i)));
            spkbinary=logical(data.linear_track.left{1,i}.dataspks(:,end));
            plot(x,ts,'.k');
            axis tight
            hold on;box off; axis off
            scatter(x(spkbinary),ts(spkbinary),20,...
                interp1(data.lfp.ts,...
                data.lfp.theta_phase(trodeID,:),ts(spkbinary)','linear'),'filled');
            colormap(ax,c)
            
            % PLOT EACH LAPS RATEMAP WITH THE OVERALL FIRING RATE
            % SUPERIMPOSED
            p(5, ns).select();
            ax = gca;
            laps=reshape([ data.linear_track.left{1,i}.maps{:}],...
                [],length(data.linear_track.left{1,i}.maps));
            imagesc(flipud(imrotate(laps,90)));hold on;set(gca,'Ydir','Normal')
            ratemap= data.ratemap(ns,:,i);
            plot(rescale(ratemap,1,size(laps,2)),'LineWidth',2, 'color','w');
            set(gca,'XTickLabel',[]);
            axis tight
            hold on;box off;
            colormap(ax,jet(255))
            title(sprintf('Info Content: %4.2f Peak Rate: %4.2f',...
                data.results.InformationContent(ns,i),data.results.PeakRate(ns,i)))
            
            % PHASE BY POSITION
            p(6, ns).select();
            ax = gca;
            x=rescale(x,0,1);
            phase=interp1(data.lfp.ts,...
                data.lfp.theta_phase(trodeID,:),ts(spkbinary)','linear');
            scatter([x(spkbinary);x(spkbinary)],[phase';phase'+2*pi]*180/pi,15,'Filled','k');
            axis tight;box off;axis off
            xlim([min(x) max(x)])

            for f=1:length(data.linear_track.left{1,i}.fields)
                rho(f)=data.linear_track.left{1,i}.fields{f}.ThPrecess.circLinCorr;
                pval(f)=data.linear_track.left{1,i}.fields{f}.ThPrecess.pval;
            end
            a3=[rho;pval];
            title(sprintf(repmat(' rho:%4.2f p:%4.2f,',1,length(rho)),a3(:)'))
            clear rho pval
            
            %                 hold on
            %                 % plot correlation lines
            %                 plot_corr_line(x(spkbinary),length(ratemap),data.linear_track.left{1,i})
            
            p(7, ns).select();
            ax = gca;
            ts=data.linear_track.left{1,i}.dataspks(data.linear_track.left{1,i}.dataspks(:,end)==0,1);
            x=data.linear_track.left{1,i}.dataspks(data.linear_track.left{1,i}.dataspks(:,end)==0,2);
            spkts=data.linear_track.left{1,i}.dataspks(data.linear_track.left{1,i}.dataspks(:,end)==1,1);
            [~,I]=unique(ts);
            ts=ts(I);
            x=x(I);
            Phasemap(ts,x,spkts,data.lfp,length(ratemap),trodeID,data.samplerate)
            colormap(ax,jet(255))
            
            
            % PLOT AUTOCORRELATION
            p(8, ns).select();
            ax = gca;
            y=data.autocor(ns,:,i);
            plot(y,'LineWidth',2, 'color','k');
            axis tight
            hold on;box off; axis off
            title(sprintf('Theta modulation: %4.2f',data.results.thetaindex(ns,i)))
            
           
        end
    end
end
end
% Subfunctions
function Phasemap(ts,x,spkts,lfp,bins,trodeID,samplerate)

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
h = 1.5;
h=4;
myfilter = fspecial('gaussian',[4 12]*h, h);
phasemap = imfilter([phasemap,phasemap,phasemap],myfilter,'replicate');
phasemap=phasemap(:,bins:(bins-1)*2);

pcolor(flipud(rot90(phasemap)));shading flat;box off;axis off;axis tight;
end
