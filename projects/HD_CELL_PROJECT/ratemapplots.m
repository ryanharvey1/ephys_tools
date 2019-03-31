% ratemapplots
% Plot and Save ratemaps with border score,bordermod,egomod,gridscore,info, autocor with theta index

load('D:\Projects\Multi_Region_HD\HDdata_fieldinfo.mat')

areas=fieldnames(HDdata);

for i=1:length(areas)
    maps=HDdata.(areas{i}).ratemap;
    
    for ncell=1:length(maps)
        ratemap=HDdata.(areas{i}).ratemap{ncell};
        frames=HDdata.(areas{i}).frames{ncell};  
        
        fig=figure;fig.Color=[1 1 1];
        subplot(1,2,1)
        x=median([frames(:,2),frames(:,4)],2);
        y=median([frames(:,3),frames(:,5)],2);
        plot(x(frames(:,6)==0),y(frames(:,6)==0),'LineWidth',2,'color','k');
        hold on; axis off
        scatter(x(frames(:,6)==1),y(frames(:,6)==1),20,'r','filled');
        box off; axis image

        subplot(1,2,2)
        imAlpha=ones(size(ratemap));
        imAlpha(isnan(ratemap))=0;
        imagesc(ratemap,'AlphaData',imAlpha);
        axis xy; axis off; hold on; box off; axis image;
        colormap(gca,jet(255))
        colorbar
        %         str=sprintf(['BorderScore: ',num2str(round(HDdata.(areas{i}).borderscore(ncell),2)),...
        %             ' BorderMod: ',num2str(round(max(HDdata.(areas{i}).bordermodulation(ncell,:)),2)),'\n',...
        %             ' EgoMod: ',num2str(round(max(HDdata.(areas{i}).egocentricmodulation(ncell,:)),2)),...
        %             ' GridScore: ',num2str(round(HDdata.(areas{i}).gridscore(ncell),2)),...
        %             ' InfoContent: ',num2str(round(HDdata.(areas{i}).informationContent(ncell),2))]);
        str=sprintf([' nfields: ',num2str(round(max(HDdata.(areas{i}).nfields(ncell,:)),2)),...
            ' GridScore: ',num2str(round(HDdata.(areas{i}).gridscore(ncell),2)),...
            ' InfoContent: ',num2str(round(HDdata.(areas{i}).informationContent(ncell),2))]);
        suptitle(str)
         
        id=strsplit(HDdata.(areas{i}).id{ncell},'.');
        id=id{1};
        
        if ~exist(['D:\Projects\Multi_Region_HD\ratemaps',filesep,areas{i}],'dir')
            mkdir(['D:\Projects\Multi_Region_HD\ratemaps',filesep,areas{i}])
        end
        
        print(fig,'-dpng', '-r150',['D:\Projects\Multi_Region_HD\ratemaps',filesep,areas{i},filesep,id,'.png'])
        close all
    end
end