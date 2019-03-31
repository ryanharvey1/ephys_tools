cd D:\ClarkP30_Recordings\CFCResults\CFCResults

control={'RH13','RH14','LS21','LS23','LE2821','LE2823','LEM3116','LEM3120'};
pae={'RH11','RH16','LS17','LS19','LE2813'};

rats=dir('*.mat');
rats={rats.name};

data1=load(rats{1});
paeamps=data1.cfc{1, 1}.LGtortCFC.MeanAmp ./sum(data1.cfc{1, 1}.LGtortCFC.MeanAmp);

data2=load(rats{2});
WTamps=data2.cfc{1, 1}.LGtortCFC.MeanAmp./ sum(data2.cfc{1, 1}.LGtortCFC.MeanAmp);

subplot(1,2,2)
bar([paeamps paeamps])
subplot(1,2,1)
bar([WTamps WTamps])

set(ax,'FontSize',20,'FontWeight','bold','LineWidth',2,'box','off')

cd d:\Users\BClarkLab\Desktop\Laura Temp\SFN2018
print(gcf,'-dmeta',['d:\Users\BClarkLab\Desktop\Laura Temp\SFN2018',...
    filesep,'SPECMI_FIG.emf'])

%PRINT SPECT MI

 imagesc(thetarange,gammarange,CFC.modindex);
    axis xy
    
    title(['Data MI  (maxmod = ' num2str(max(CFC.modindex(:)))])
    xlabel('Theta Frequency (Hz)'), ylabel('Gamma Frequency (Hz)')
    
%     subplot(1,2,2)
%     imagesc(thetarange,gammarange,CFC.modindex_shuffled);
%     axis xy
%     
%     title(['Surrogate MI  (maxmod = ' num2str(max(CFC.modindex_shuffled(:)))])
%     xlabel('Theta Frequency'), ylabel('Gamma Frequency')
     
     
    tempMI=data.cfc{1, 1}.LGtortCFC.modindex; 
    tempThres=data.cfc{1, 1}.LGtortCFC.MIthres;
    
    if tempMI<tempThres
        continue
    end
    
    if contains(rats(i),control) && .0016+.0002 >= tempMI &&  .0016-.0002 <= tempMI %thresholds obtained from session analysis (mean+-std)
        disp(['Control Sess to check ',rats{i}])
    elseif contains(rats(i),pae) && .0008+.0001 >= tempMI && .0008-.0001<=tempMI %thresholds obtained from session analysis (mean+-std)
        disp(['PAE Sess to check ',rats{i}])
    end

