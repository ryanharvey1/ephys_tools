%CumulativeFreq 
% Messy code used to create some figures for RSA 2017
close all
measure=1;
figure;
data1=controlCD(:,1);
data2=PAECD(:,1);
% data1(data1==0)=[]; data2(data2==0)=[];
[p,h,stats] = ranksum(data1,data2)

[f,x] = ecdf(data1);
p4=plot(x,f); x1=x;
set(p4,'LineWidth',4,'Color','k')
hold on
[f,x] = ecdf(data2);
p4_2=plot(x,f); x2=x;
set(p4_2,'LineWidth',4,'Color','r')
box off
legend([p4 p4_2],'Control','PAE','FontSize',12,'Location','best')
xlabel(AllVariableNames(measure));ylabel('Cumulative Frequency')

[h,p,k]=kstest2(x1,x2)
[p,h,stats] = ranksum(x1,x2)
%%
% [p,h,stats] = signrank(x1,x2)
% [p,h,stats] = signtest(x1,x2)

figure(5); subplot(1,2,1); h1=histogram(data1,10);
xlim([0 22]); box off;
title('Control'); ylabel('Number of Fields')
subplot(1,2,2); h2=histogram(data2,10);
xlim([0 22]); box off;
title('PAE');ylabel('Number of Fields');
set(h1,'FaceColor',[0.247 0.247 0.247])
set(h2,'FaceColor',[0.247 0.247 0.247])


Controldisp=histcounts(data1,22);
PAEdisp=histcounts(data2,22);

% normalize from 0-1 WORKING*********!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Controldisp= (Controldisp - min(Controldisp)) / ( max(Controldisp) - min(Controldisp) );
  PAEdisp= (PAEdisp - min(PAEdisp)) / ( max(PAEdisp) - min(PAEdisp) );
    
h1=area(Controldisp); hold on; h2=area(PAEdisp);
h1=bar(Controldisp); hold on; h2=bar(PAEdisp);

% Controldisp/sum(Controldisp)
% PAEdisp/sum(PAEdisp)



figure(5); subplot(1,2,1); h1=histogram(RowCorrelations_control,10);
% xlim([0 22]);
box off;
title('Control'); ylabel('Number of Fields')
subplot(1,2,2); h2=histogram(RowCorrelations_PAE,10);
% xlim([0 22]);
box off;
title('PAE');ylabel('Number of Fields');
set(h1,'FaceColor',[0.247 0.247 0.247])
set(h2,'FaceColor',[0.247 0.247 0.247])

[p,h,stats] = ranksum(RowCorrelations_control,RowCorrelations_PAE)

[h,p,ci,stats] = ttest2(diagCorr_control,diagCorr_PAE)
[p,h,stats] = ranksum(diagCorr_control,diagCorr_PAE)

    % Mean Half-width
    half1=decordistALLnorm(1,:); half1=half1';
    half2=decordistALLnorm(2,:); half2=half2';
%     halfwidth1=median(half1(half1>0));
%     halfwidth2=median(half2(half2>0));
halfwidth1=.5;
halfwidth2=.5;
    
    x=linspace(-120,120,43);
    x1=interp1(half1(1:22,1),x(1,1:22),halfwidth1,'linear');
    x2=interp1(half1(22:end,1),x(1,22:end),halfwidth1,'linear');
    
    x11=interp1(half2(1:22,1),x(1,1:22),halfwidth2,'linear');
    x22=interp1(half2(22:end,1),x(1,22:end),halfwidth2,'linear');

    figure(6);
    p1=plot(linspace(-120,120,43),decordistALLnorm(1,:),'color','k');
    hold on
%     plot([x1,x2],[halfwidth1;halfwidth1],'color','k')
    hold on
    p2=plot(linspace(-120,120,43),decordistALLnorm(2,:),'color','r');
    hold on
%     plot([x11,x22],[halfwidth2;halfwidth2],'color','r')
    xlabel('Distance (cm)')
    ylabel('Correlation')
    title('Decorrelation Distance')
    legend([p1 p2],'Control','PAE','FontSize',12)
    
    controlHalfWidth=x2-x1;
    PAEHalfWidth=x22-x11;

