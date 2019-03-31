% HD_Project_TempAutoCorr

clear;clc;close all
addpath('/Users/RyanHarvey/GoogleDrive/MatlabDir/BClarkToolbox/Analysis')

path='/Volumes/SAMSUNGUSB/Ben_HDCProject';
cd(path)
folders=dir;
folders={folders.name}';
folders(contains(folders,'._'))=[];
folders(contains(folders,'.'))=[];
folders(contains(folders,'SpikePositionFiles'))=[];


for a=1:length(folders)
    cd(['/Volumes/SAMSUNGUSB/Ben_HDCProject/',folders{a}])
    
    filenames=dir('*.ts_R');
    filenames={filenames.name}';
    filenames(contains(filenames,'._'))=[];
    
    for i=1:length(filenames)
        % OPEN TEXT FILE AND CLEAN DATA
        disp(['Running:  ',filenames{i}])
        fileID = fopen(filenames{i},'r');
        dataArray = textscan(fileID,'%f%*s%[^\n\r]','Delimiter',' ','MultipleDelimsAsOne',...
            true,'TextType','string','EmptyValue',NaN,'HeaderLines',2-1,'ReturnOnError',false,'EndOfLine','\r\n');
        fclose(fileID);
        ts=[dataArray{1:end-1}];
        clear dataArray fileID
        
        %% FIND NEGATIVE NUMBERS
        nNeg(i,:)=sum(ts<0);
        ts(ts<0)=[];
        %% CALC THETA MODULATION
%         ts=ts*1;
        [thetaindex(i),peak(i),cor(i,:),lag] = thetamodulation(ts/10000000);
        
        %    figure;plot(([-100:-1,0,1:100]*5), cor(i,:),'k'); hold on
        %        set(gca,'YTickLabel',[],'YTick',[],'XMinorTick','on','YMinorTick','off','LineWidth',1)
        %        line([0 0], ylim, 'linestyle', ':', 'color', [.7 .7 .7]);
        %        set(gca, 'fontsize', 20, 'box', 'off');
        %        title(['Theta Index= ',num2str(thetaindex(i)),' Freq= ',num2str(peak(i))])
    end
    %% SORT DATA
    [thetaindex,I]=sort(thetaindex,'ascend');
    cor=cor(I,:);
    
    %% Save into data
    
    currentarea=folders{a};
    index = find(isletter(currentarea), 1);
    currentarea  = currentarea(index:end);
    
    currentarea=strsplit(currentarea,'_');
    data.(currentarea{1}).autocorr=cor;
    data.(currentarea{1}).thetaindex=thetaindex;
    
end

%% GRAPH autocorrs sorted by thetaindex arranged by brain area

folders=fieldnames(data);

autocor_fig=figure; autocor_fig.Color=[1 1 1];

for a=1:length(folders)
    subplot(2,2,a)
    imagesc(data.(folders{a}).autocorr);
    colormap jet; box off;
    ylabel('Cells')
    set(gca,'XMinorTick','on','YMinorTick','off',...
        'LineWidth',1,'box','off','XTick',linspace(5,201,5),'XTickLabel',...
        [linspace(-500,-1,2),0,linspace(1,500,2)],'FontSize',20,'FontWeight','bold')
    title(folders{a})
end

%% SAVE
print(autocor_fig,'-dpng', '-r600',['/Volumes/SAMSUNGUSB',filesep,'autocor_fig.png'])

%% GRAPH All tuning curves in one matrix sorted by RLength

ALL_ThetaIdx=[];
ALL_autocorr=[zeros(1,size(data.(folders{1}).autocorr,2))];
for a=1:length(folders)
    ALL_ThetaIdx=[ALL_ThetaIdx;data.(folders{a}).thetaindex'];
    ALL_autocorr=[ALL_autocorr;data.(folders{a}).autocorr];
end
ALL_autocorr(1,:)=[];

[SortedThetaIdx,I]=sort(ALL_ThetaIdx);
Sortedautocor=ALL_autocorr(I,:);

ALLautocor_Fig=figure; ALLautocor_Fig.Color=[1 1 1];

imagesc(Sortedautocor);
colormap jet; box off;
ylabel('Cells')
set(gca,'XMinorTick','on','YMinorTick','off',...
    'LineWidth',1,'box','off','XTick',linspace(5,201,5),'XTickLabel',...
    [linspace(-500,-1,2),0,linspace(1,500,2)],'FontSize',20,'FontWeight','bold')
title('ATN,MEC,PaS,PoS')
%% SAVE
print(ALLautocor_Fig,'-dpng', '-r600',['/Volumes/SAMSUNGUSB',filesep,'ALLautocor_Fig.png'])

%% SAVE TO CSV
csvwrite('/Volumes/SAMSUNGUSB/Sortedautocor.csv',Sortedautocor)
