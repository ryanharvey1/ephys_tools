% Theta_modulation
clear;clc;close all
addpath ('/Users/ryanharvey/GoogleDrive/MatlabDir/CircStat2012a','/Users/ryanharvey/GoogleDrive/MatlabDir/BClarkToolbox/Analysis')
% load('/Users/ryanharvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/NewOCC_Map_workspace2.mat')
% load('/Users/ryanharvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/data4_SPLITUPBYREGION');
load('/Users/RyanHarvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/NewPipelineResults.mat')
FigureLocation='/Users/RyanHarvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/Figures_For_Paper';
%%

[thetaindex,cor,peak,modulated]=computethetamod(results.control.data);
thetamod.control.thetaindex=thetaindex;
thetamod.control.cor=cor;
thetamod.control.peak=peak;
thetamod.control.shuffmod=modulated;

Group1=[thetaindex',peak'];

[s,I]=sort(thetaindex);
cor=cor(I,:);
figure;imagesc(cor);colormap jet;hold on;box off;axis off;
plot(rescale(-sum(cor),1,size(cor,1)),'w','LineWidth',3)

[thetaindex,cor,peak,modulated]=computethetamod(results.tilted.data );
thetamod.tilted.thetaindex=thetaindex;
thetamod.tilted.cor=cor;
thetamod.tilted.peak=peak;
thetamod.tilted.shuffmod=modulated;

Group2=[thetaindex',peak'];

[s,I]=sort(thetaindex);
cor=cor(I,:);
figure;imagesc(cor);colormap jet;hold on;box off;axis off;
plot(rescale(-sum(cor),1,size(cor,1)),'w','LineWidth',3)

[ AllStats ] = CDFplots(Group1(:,1),Group2(:,1),{'Control','Tilted'},{'Theta Modulation'},2)

%% plot cells that passed shuffled distribution
thetamod.control.thetaindex(logical(thetamod.control.shuffmod)) 
cor=thetamod.control.cor(logical(thetamod.control.shuffmod),:);

[s,I]=sort(thetamod.control.thetaindex(logical(thetamod.control.shuffmod)));
cor=cor(I,:);
figure;imagesc(cor);colormap jet;hold on;box off;axis off;
plot(rescale(-sum(cor),1,size(cor,1)),'w','LineWidth',3)

thetamod.tilted.thetaindex(logical(thetamod.tilted.shuffmod)) 
cor=thetamod.tilted.cor(logical(thetamod.tilted.shuffmod),:);

[s,I]=sort(thetamod.tilted.thetaindex(logical(thetamod.tilted.shuffmod)));
cor=cor(I,:);
figure;imagesc(cor);colormap jet;hold on;box off;axis off;
plot(rescale(-sum(cor),1,size(cor,1)),'w','LineWidth',3)

% compare peak theta freq
[ AllStats ] = CDFplots(thetamod.control.peak(logical(thetamod.control.shuffmod))',thetamod.tilted.peak(logical(thetamod.control.shuffmod))',{'Control','Tilted'},{'Theta Frequency'},2)

print(figure(3), '-dpng', '-r600',['/Users/RyanHarvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/Figures_For_Paper',filesep,'Theta_Freq.png'])
print(figure(4), '-dpng', '-r600',['/Users/RyanHarvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/Figures_For_Paper',filesep,'ThetaModulation.png'])



%%

% stuct4shuff=results.tilted.data;
% 
%  f = fieldnames(results.control.data);
%  for i = 1:length(f)
%     stuct4shuff.(f{i}) = results.control.data.(f{i});
%  end
%  
% 
% 
%  ids=fieldnames(stuct4shuff);
%  j=1;
%  for i=1:5:length(ids)
%      tempframes=stuct4shuff.(ids{i});
%      for x=1:400
%          tempframes=circshift(tempframes);
%          spk=tempframes(tempframes(:,5)==1,1);
%          [thetaindex(j),peak(j),cor(j,:),lag] = thetamodulation(spk/1000000);
%          j=j+1;
%      end
%      
%  end
% 




% print(figure(4), '-dpng', '-r600',['/Users/RyanHarvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/Figures_For_Paper',filesep,'ThetaModulation_Freq.png'])

% %% splitupbyregion
% reidx.sessXregionC=sessXregionC;
% reidx.sessXregionT=sessXregionT;
% 
% DispVarNames={'PeakRate','nSpikes','OverallFR','NumbActiveBins','sparsity','InformationContent','Coherence','Field2Wall','borderScore','FieldWidth','infieldFR','outfieldFR','E','c','p','Theta Modulation','Theta Frequency'};
% 
% resultsC(:,:,[2:5])=[];
% resultsT(:,:,[2:5])=[];
% 
% resultsC=[resultsC,Group1];
% resultsT=[resultsT,Group2];
% 
% [ResultsC]=analyze(sessionC,sessXregionC,resultsC,reidx);
% [ResultsT]=analyze(sessionT,sessXregionT,resultsT,reidx);
% ca1c=ResultsC.ca1;
% ca3c=ResultsC.ca3;
% ca1t=ResultsT.ca1;
% ca3t=ResultsT.ca3;
% 
% % [ CvC ] = ScatterBox(ca1c(:,:,1),ca3c(:,:,1),{'Con CA1','Con CA3'},DispVarNames,1)
% % [ TvT ] = ScatterBox(ca1t(:,:,1),ca3t(:,:,1),{'Tilt CA1','Tilt CA3'},DispVarNames,1)
% % 
% % [ CvTca1 ] = ScatterBox(ca1c(:,:,1),ca1t(:,:,1),{'Con CA1','Tilt CA1'},DispVarNames,1)
% % [ CvTca3 ] = ScatterBox(ca3c(:,:,1),ca3t(:,:,1),{'Con CA3','Tilt CA3'},DispVarNames,1)
% % 
% % [ CvT ] = ScatterBox([ca1c(:,:,1);ca3c(:,:,1)],[ca1t(:,:,1);ca3t(:,:,1)],{'Con','Tilt'},DispVarNames,1)
% % 
% 
% 
% [ CvC ] = ScatterBox(ca1c(:,[16:17],1),ca3c(:,[16:17],1),{'Con CA1','Con CA3'},DispVarNames(16:17),2)
% [ TvT ] = ScatterBox(ca1t(:,[16:17],1),ca3t(:,[16:17],1),{'Tilt CA1','Tilt CA3'},DispVarNames(16:17),2)
% 
% 
% 
% function [results]=analyze(sess,sessXregion,data,reidx)
% 
% sessXregion=sessXregion(~cellfun('isempty',sessXregion(:,1)),:);
% sessXregion=sessXregion(~cellfun('isempty',sessXregion(:,2)),:);
% 
% sess=sess(~cellfun('isempty',sess(:,1)),:);
% 
% sessXregion(ismember(sessXregion(:,2),'dentate'),2)={'CA3'};
% 
% CA1=sessXregion(ismember(sessXregion(:,2),'CA1'),1);
% CA3=sessXregion(ismember(sessXregion(:,2),'CA3'),1);
% 
% idxca1=ismember(sess(:),CA1);
% idxca3=ismember(sess(:),CA3);
% 
% results.ca1=data(idxca1,:,:);
% results.ca3=data(idxca3,:,:);
% end
%%
function [thetaindex,cor,peak,modulated]=computethetamod(data)
ids=fieldnames(data);
j=1;
for i=1:5:length(ids)
    tempframes=data.(ids{i});
    spk=tempframes(tempframes(:,5)==1,1);
    [thetaindex(j),peak(j),cor(j,:),lag] = thetamodulation(spk/1000000);
    %        cortemp=smooth(cortemp,5)';
    %        cor(j,:)=cortemp;
    %    figure;plot(([-50:-1,0,1:50]*10), cor,'k'); hold on
    %    set(gca,'YTickLabel',[],'YTick',[],'XMinorTick','on','YMinorTick','off','LineWidth',1)
    %    line([0 0], ylim, 'linestyle', ':', 'color', [.7 .7 .7]);
    %    set(gca, 'fontsize', 20, 'box', 'off');
    %    title(['Theta Index= ',num2str(thetaindex(j)),' Freq= ',num2str(peak)])
    %
    
    % Pass Shuff?
    for ishuff=1:400
%         tempframes(:,5)=circshift(tempframes(:,5),randi(length(tempframes)));
        tempframes(:,5)=tempframes(randperm(length(tempframes)),5);

        spk=tempframes(tempframes(:,5)==1,1);
        [thetaindexshuff(ishuff,1),~,~,~] = thetamodulation(spk/1000000);
    end

    if thetaindex(j)>=prctile(thetaindexshuff,95)
        modulated(j,1)=1;
    else
        modulated(j,1)=0;
    end
    clear thetaindexshuff spk tempframes
    
    j=j+1;
end
end