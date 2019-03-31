% FiringbyCue
% Comment 1: It appears that the place fields in tilted mice cluster near
% the cue card (and not just any border). It would be helpful to quantify
% this apparent clustering near the cue card. Such figures could replace
% of some of the B-E figures, which are redundant and thus could be
% covered in the text but not a figure.
% QUANTIFY LOCATIONS OF FIELDS TO SEE PORPORTION AROUND CUE 

clear;clc;close all
addpath('/Users/RyanHarvey/GoogleDrive/MatlabDir/BClarkToolbox/Analysis')
load('/Users/ryanharvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/NewOCC_Map_workspace2.mat')
FigureLocation='/Users/RyanHarvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/Figures_For_Paper';
[resultsC]=analyze(ratemapC,1);
[resultsT]=analyze(ratemapT,1);

for i=1:4
    if i==1; disp('T');end; if i==2; disp('R'); end; if i==3; disp('B'); end; if i==4; disp('L');end;
    [h,p, chi2stat,df] = prop_test([resultsC(i),resultsT(i)] , [length(fieldnames(ratemapC))/5 length(fieldnames(ratemapT))/5], 0)
end
% [table,chi2,p] = crosstab(resultsC,resultsT)

       % Observed data
       n1 = 51; N1 = 8193;
       n2 = 74; N2 = 8201;
       x1 = [repmat('a',N1,1); repmat('b',N2,1)];
       x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
       [tbl,chi2stat,pval] = crosstab(x1,x2)
       
       
       n1 = round(resultsC(1)*100); N1 = 100;
       n2 = round(resultsC(2)*100); N2 = 100;
       n3 = round(resultsC(3)*100); N3 = 100;
       n4 = round(resultsC(4)*100); N4 = 100;
       x1 = [repmat('a',N1,1); repmat('b',N2,1); repmat('c',N3,1); repmat('d',N4,1)];
       x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1); repmat(1,n3,1); repmat(2,N3-n3,1); repmat(1,n4,1); repmat(2,N4-n4,1)];
       [tbl,chi2stat,pval] = crosstab(x1,x2)
       
       n1 = round(resultsT(1)*100); N1 = 100;
       n2 = round(resultsT(2)*100); N2 = 100;
       n3 = round(resultsT(3)*100); N3 = 100;
       n4 = round(resultsT(4)*100); N4 = 100;
       x1 = [repmat('a',N1,1); repmat('b',N2,1); repmat('c',N3,1); repmat('d',N4,1)];
       x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1); repmat(1,n3,1); repmat(2,N3-n3,1); repmat(1,n4,1); repmat(2,N4-n4,1)];
       [tbl,chi2stat,pval] = crosstab(x1,x2)

       n1 = round(resultsC(1)*100); N1 = 100;
       n2 = round(resultsC(2)*100); N2 = 100;
       n3 = round(resultsC(3)*100); N3 = 100;
       n4 = round(resultsC(4)*100); N4 = 100;
       n5 = round(resultsT(1)*100); N5 = 100;
       n6 = round(resultsT(2)*100); N6 = 100;
       n7 = round(resultsT(3)*100); N7 = 100;
       n8 = round(resultsT(4)*100); N8 = 100;
       x1 = [repmat('a',N1,1); repmat('b',N2,1); repmat('c',N3,1); repmat('d',N4,1); repmat('e',N5,1); repmat('f',N6,1); repmat('g',N7,1); repmat('h',N8,1)];
       x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1); repmat(1,n3,1); repmat(2,N3-n3,1); repmat(1,n4,1); repmat(2,N4-n4,1);...
           repmat(1,n5,1); repmat(2,N5-n5,1); repmat(1,n6,1); repmat(2,N6-n6,1); repmat(1,n7,1); repmat(2,N7-n7,1); repmat(1,n8,1); repmat(2,N8-n8,1)];
       [tbl,chi2stat,pval,label] = crosstab(x1,x2)
       

%         print(figure(1),'-dpng', '-r600',[FigureLocation,filesep,'controlFieldLocation2.png'])
%         print(figure(2),'-dpng', '-r600',[FigureLocation,filesep,'tiltedFieldLocation2.png'])

function [results]=analyze(ratemap,normfield)
mapC=[];
sessions=fieldnames(ratemap);
sessions=sessions([1:5:length(sessions)],1);
results=zeros(length(sessions),4);
for i=1:length(sessions)
    % pull in matrix
    mapC=ratemap.(sessions{i});
    % find field
    [field,~]=FindFF2D(mapC);
    if normfield==0
        ALL_fields(:,:,i)=field;
    elseif normfield==1
        % isolate field
        tempmap=mapC;
        tempmap(~field)=0;
        % normalize field 0 to 1
        tempmap=tempmap-min(tempmap(:));
        ALL_fields(:,:,i)=tempmap./max(tempmap(:));
    end
    
    T=sum(sum(flipud(fliplr(triu(fliplr(triu(flipud(field))))))))/sum(sum(field));
    R=sum(sum(flipud(triu(flipud(triu(field))))))/sum(sum(field));
    B=sum(sum(fliplr(triu(fliplr(triu(field))))))/sum(sum(field));
    L=sum(sum(fliplr(flipud(triu(flipud(triu(fliplr(field))))))))/sum(sum(field));
    
    %     figure;imagesc(fliplr(triu(fliplr(triu(field)))));axis xy
    %     figure;imagesc(flipud(triu(flipud(triu(field)))));axis xy
    %     figure;imagesc(flipud(fliplr(triu(fliplr(triu(flipud(field)))))));axis xy
    %     figure;imagesc(fliplr(flipud(triu(flipud(triu(fliplr(field)))))));axis xy
    %     figure;imagesc(ALL_fields(:,:,i));axis xy
    %     close all
    
    [M,I]=max([T,R,B,L]);
    results(i,I)=1;
end
results=sum(results)/length(results);

ALLfields=sum(ALL_fields,3);
upsamRateMap=PerfectCircRateMap(ALLfields,0);
imAlpha=ones(size(upsamRateMap));
imAlpha(isnan(upsamRateMap))=0;
figC=figure;imagesc(upsamRateMap,'AlphaData',imAlpha);
figC.Color=[1 1 1];
box off
axis image
axis off
axis xy
colormap jet
% shading interp

% hold on;  plot(1:1000,1:1000,'w','LineWidth',6);plot(1000:-1:1,1:1000,'w','LineWidth',6)
% txt=text(500,250,num2str(controlresults(1)),'FontSize',20,'Color',[1 1 1]);
% txt=text(750,500,num2str(controlresults(2)),'FontSize',20,'Color',[1 1 1]);
% txt=text(500,750,num2str(controlresults(3)),'FontSize',20,'Color',[1 1 1]);
% txt=text(250,500,num2str(controlresults(4)),'FontSize',20,'Color',[1 1 1]);

end

% clear tempmap ALL_fields results


%
% mapC=[];
% sessions=fieldnames(ratemapT);
% sessions=sessions([1:5:495],1);
% results=zeros(length(sessions),4);
% for i=1:length(sessions)
%     % pull in matrix
%     mapC=ratemapT.(sessions{i});
%     % find field
%     [field,fieldwidth]=FindFF2D(mapC);
%     % isolate field
%     tempmap=mapC;
%     tempmap(~field)=0;
%     ALL_fields(:,:,i)=tempmap;
%     % normalize field 0 to 1
% %     tempmap=tempmap-min(tempmap(:));
% %     ALL_fields(:,:,i)=tempmap./max(tempmap(:));
%
%     T=sum(sum(fliplr(triu(fliplr(triu(field))))))/sum(sum(field));
%     R=sum(sum(flipud(triu(flipud(triu(field))))))/sum(sum(field));
%     B=sum(sum(flipud(fliplr(triu(fliplr(triu(flipud(field))))))))/sum(sum(field));
%     L=sum(sum(fliplr(flipud(triu(flipud(triu(fliplr(field))))))))/sum(sum(field));
%
%     [M,I]=max([T,R,B,L]);
%     results(i,I)=1;
% end
% tiltedresults=sum(results)./length(results)
% % [h,p, chi2stat,df] = prop_test([54 100] , [422 1204], 0);
%
% ALLfields=sum(ALL_fields,3);
% upsamRateMap=PerfectCircRateMap(ALLfields,0);
% imAlpha=ones(size(upsamRateMap));
% imAlpha(isnan(upsamRateMap))=0;
% figT=figure;imagesc(upsamRateMap,'AlphaData',imAlpha);
% figT.Color=[1 1 1];
% box off
% axis image
% axis off
% axis xy
% colormap jet
% shading interp
%
% hold on;  plot(1:1000,1:1000,'w','LineWidth',2);plot(1000:-1:1,1:1000,'w','LineWidth',2)
% % txt=text(500,250,num2str(tiltedresults(1)),'FontSize',20,'Color',[1 1 1]);
% % txt=text(750,500,num2str(tiltedresults(2)),'FontSize',20,'Color',[1 1 1]);
% % txt=text(500,750,num2str(tiltedresults(3)),'FontSize',20,'Color',[1 1 1]);
% % txt=text(250,500,num2str(tiltedresults(4)),'FontSize',20,'Color',[1 1 1]);
% %



%     T=fliplr(triu(fliplr(triu(ALLfields))));
%     R=flipud(triu(flipud(triu(ALLfields))));
%     B=flipud(fliplr(triu(fliplr(triu(flipud(ALLfields))))));
%     L=fliplr(flipud(triu(flipud(triu(fliplr(ALLfields))))));
%
%     figure;imagesc(T);
% figure;imagesc(R);
% figure;imagesc(B);
% figure;imagesc(L);

% tempmap=sum(All_map,3);
%
%     % cut up the matrix pie and store them
%     T=fliplr(triu(fliplr(triu(tempmap))));
%     R=flipud(triu(flipud(triu(tempmap))));
%     B=flipud(fliplr(triu(fliplr(triu(flipud(tempmap))))));
%     L=fliplr(flipud(triu(flipud(triu(fliplr(tempmap))))));

% figure;imagesc(T);
% figure;imagesc(R);
% figure;imagesc(B);
% figure;imagesc(L);

% figure;imagesc(mapC);
% hold on; plot(bound(:,1),bound(:,2))
% hold on; plot(bound(:,1),bound(:,2),'r')
% hold on;  plot(1:25,1:25,'r')
% hold on;  plot(25:-1:1,1:25,'r')
%
% T=fliplr(triu(fliplr(triu(field))));
% figure;imagesc(T)
%
% R=flipud(triu(flipud(triu(field))));
% figure;imagesc(R)
%
% B=flipud(fliplr(triu(fliplr(triu(flipud(field))))));
% figure;imagesc(B)
%
% L=fliplr(flipud(triu(flipud(triu(fliplr(field))))));
% figure;imagesc(L)




% mapT=[];
% sessions=fieldnames(ratemapT);
% sessions=sessions([1:5:225],1);
% for i=1:length(sessions)
%     mapT=[mapT;ratemapT.(sessions{i})];
% end


