% HD_angular_linear_velocity
if ~exist('HDdata','var')
    load('D:\Projects\Multi_Region_HD\HDdata')
end
areas=fieldnames(HDdata);


bins=33;
for a=1:length(areas)
   for cells=1:length(HDdata.(areas{a}).frames) 
       
       frames=HDdata.(areas{a}).frames{cells};
       frames_w_spikes=frames;
       frames(frames(:,6)==1,:)=[];
       
       anglevel=insta_angvel(rad2deg(frames(:,10)),60);

       x=median([frames(:,2),frames(:,4)],2);
       y=median([frames(:,3),frames(:,5)],2);
       
       x=smoothdata(x,'loess',10);
       y=smoothdata(y,'loess',10);
       
       if contains(areas{a},'ATN')
           mazesize=70;
       else
           mazesize=140;
       end
       [vel_cmPerSec,vel_abs,pixelDist]=InstaVel([x,y],'no',mazesize,60);
       
       angle_edges=linspace(-300,300,bins);
       linearvel_edges=linspace(0,100,bins);
       occ=histcounts2(vel_cmPerSec,anglevel,linearvel_edges,angle_edges);
       
       lowocc=occ==0;
       
       % spike map
       vel_cmPerSec=interp1(frames(2:end,1),vel_cmPerSec,frames_w_spikes(frames_w_spikes(:,6)==1,1));
       anglevel=interp1(frames(2:end,1),anglevel,frames_w_spikes(frames_w_spikes(:,6)==1,1));

       spikemap=histcounts2(vel_cmPerSec,anglevel,linearvel_edges,angle_edges);
       
       AHV_LV_map=spikemap./(occ*60);
       
       AHV_LV_map(isnan(AHV_LV_map) | isinf(AHV_LV_map))=0;
       AHV_LV_map(lowocc)=NaN;
       
       filtWidth = [5 5]; filtSigma = 1;
       imageFilter=fspecial('gaussian',filtWidth,filtSigma);
       AHV_LV_map = nanconv(AHV_LV_map,imageFilter, 'nanout');
       
       fig=figure;fig.Color=[1 1 1];
       pcolor(AHV_LV_map)
       colormap jet
       box off
       shading flat
       colorbar;
       set(gca,'XTick',linspace(1,size(AHV_LV_map,2),2),...
           'XTickLabel',[min(angle_edges) max(angle_edges)],...
           'YTick',linspace(1,size(AHV_LV_map,1),2),...
           'YTickLabel',linspace(min(linearvel_edges),max(linearvel_edges),2),...
           'FontSize',20,'FontWeight','bold','LineWidth',3)
        xlabel('Angular Velocity (deg/sec)')
        ylabel('Linear Velocity (cm/sec)')
        
        id=HDdata.(areas{a}).id{cells};
        disp(id)
        print(fig,'-dpng', '-r100',['C:\Users\ryanh\Dropbox\school work\UNM\Lab\Projects\Multi_Region_HDcells\Figures\AHV_LV',filesep,areas{a},'_',id,'.png'])
        close      
   end
    
end