anal_dir = '~/Dropbox/Analysis_Thalamus/';
figDir = [anal_dir 'Figures/AdnPostsub/'];

bin = 0.5;
nbBins = 100;

xcWAll = [];
xcSAll = [];
xcRAll = [];

cellTypeAll = [];

cch = [];
mwf = [];
spkW = [];
dayIx = [];
shankIx = [];
neuIx = [];
cIx = [];
poHd = [];

for ii=1:length(datasets)

    
    [dumy fbasename dumy] = fileparts(datasets{ii});
    disp(fbasename)
    load([datasets{ii} '/Analysis/SpikeData.mat']);
    load([datasets{ii}  '/Analysis/HDCells.mat']);
    load([datasets{ii}  '/Analysis/BehavEpochs.mat']);
    load([datasets{ii}  '/Analysis/GeneralInfo.mat']);
    load([datasets{ii}  '/Analysis/SpikeWaveF.mat']);
    
    hdC = hdCellStats(:,end);

    remEp = LoadEpoch([datasets{ii} '/' fbasename],'REM');
    swsEp = LoadEpoch([datasets{ii} '/' fbasename],'SWS');

    shIxPo = find(ismember(shank,shankStructure{'postsub'}));
    shIxTh = find(ismember(shank,shankStructure{'thalamus'}));
    shIxPo = shIxPo(:)';
    shIxTh = shIxTh(:)';
    
    for x=shIxTh
        rgx = Range(S{x});
        rx = length(Range(S{x}));
        for y=shIxPo
            rgy = Range(S{y});
            ry = rate(S{x});
            ryw = rate(S{y},wakeEp);
            ryr = rate(S{y},remEp);
            rys = rate(S{y},swsEp);
%             rxw = rate(S{x},wakeEp);
%             rxr = rate(S{x},remEp);
%             rxs = rate(S{x},swsEp);
            
            [h,b] = CrossCorr(rgx,rgy,bin,nbBins);
    
             rgwx = Range(Restrict(S{x},wakeEp));
             rgwy = Range(Restrict(S{y},wakeEp));
             rgrx = Range(Restrict(S{x},remEp));
             rgry = Range(Restrict(S{y},remEp));
             rgsx = Range(Restrict(S{x},swsEp));
             rgsy = Range(Restrict(S{y},swsEp));

             [hw,b] = CrossCorr(rgwx,rgwy,bin,nbBins);
             [hs,b] = CrossCorr(rgsx,rgsy,bin,nbBins);
             [hr,b] = CrossCorr(rgrx,rgry,bin,nbBins);
             [ha,b] = CrossCorr(rgx,rgy,bin,nbBins);

             xcWAll = [xcWAll hw/(ryw)];
             xcSAll = [xcSAll hs/(rys)];
             xcRAll = [xcRAll hr/(ryr)];
             xcAll = [xcAll ha/(ry)];
             cch = [cch h*rx*0.001];

            cIx = [cIx;[x y]];
            cellTypeAll = [cellTypeAll;cellType(y)];
            mwf = [mwf;meanWaveF{y}(maxIx(y),:)];
            spkW = [spkW;spkWidth(y)];
            dayIx = [dayIx;ii];
            shankIx = [shankIx;shank(y)];
            neuIx = [neuIx;cellIx(y)];
            poHd = [poHd;hdC(y)];
          
        end
    end
end

bix = b>-8 & b<8;
[ PVALS, PRED, QVALS ] = cch_conv(round(cch),21);
PVALS = PVALS(bix,:);
QVALS = QVALS(bix,:);

gbUpper = poissinv( 1 - 0.001, max( PRED(bix,:), [], 1 ) );
hiBins = bsxfun( @gt, cch( bix, : ), gbUpper ) & ( cch( bix, : ) > 0 );

%figure(1),clf,imagesc(log10(PVALS+eps)),caxis([-3 0])

%dp1 = PVALS(1:end-1,:);
%dp2 = PVALS(2:end,:);
dp1 = hiBins(1:end-1,:);
dp2 = hiBins(2:end,:);

ixe = dp1 & dp2;
ixExc = find(sum(ixe)); 
ixExcBo = (sum(ixe)); 

dp1 = QVALS(1:end-1,:);
dp2 = QVALS(2:end,:);
ixi = dp1<0.001 & dp2<0.001;
ixInh = find(sum(ixi & ~ixe)); 

xcorrPk = zeros(length(ixInh),1);
cellT = zeros(length(ixInh),1);
zVal = zeros(length(ixInh),1);

bb = b(bix);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Excitatory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xcorrPk = zeros(length(ixExc),1);
cellT = zeros(length(ixExc),1);
zVal = zeros(length(ixExc),1);

for ii=1:length(ixExc)
    
    [dumy xcorrPk(ii)] = min(PVALS(:,ixExc(ii)));
    cellT(ii) = cellTypeAll(ixExc(ii));
    z = zscore(cch(:,ixExc(ii)));
    zb = z(b>-10 & b<10);
    
    zVal(ii) = zb(xcorrPk(ii));
    
    if 0
    figure(1),clf
    subplot(1,2,1)
    plot(b,z)
    subplot(1,2,2)
    plot(mwf(ixExc(ii),:))
    title([cellT(ii) bb(xcorrPk(ii)) spkW(ixExc(ii))])
    pause
    end
end

hi = hist(bb(xcorrPk(cellT==1)),bb);
hp = hist(bb(xcorrPk(cellT==2)),bb);
hu = hist(bb(xcorrPk(cellT==0)),bb);

figure(1),clf
bar(bb,[hi;hp;hu;]','stacked')
colormap  hot
legend('INT','PYR','Undet')
if 0
    epswrite([figDir 'HDCell_MonoSyn_Hist.eps'],'size','screen')
end

ix = ixExc;
ix(bb(xcorrPk)<0) = [];
cT = cellT;
cT(bb(xcorrPk)<0) = [];

ce = zscore(cch(:,ix))';
cew = zscore( xcWAll(:,ix))';
ces = zscore( xcSAll(:,ix ))';
cer = zscore( xcRAll(:,ix))';
cerw = zscore(xcRAll(:,ix)-xcWAll(:,ix))';
cesw = zscore(xcSAll(:,ix )-xcWAll(:,ix))';

axE = (1:sum(cellT==2));
colormap jet

figure(2),clf
subplot(3,2,1)
imagesc(b,axE,ce(cT==2,:)),colorbar
caxis([-2.5 6])
xlim([-12 12])
title('PostSyn PYR, All states')
ylabel('Thal-Pos Pair #')

subplot(3,2,2)
imagesc(b,axE,cew(cT==2,:)),colorbar
caxis([-2.5 6])
xlim([-12 12])
title('Wake')
subplot(3,2,3)
imagesc(b,axE,ces(cT==2,:)),colorbar
caxis([-2.5 6])
xlim([-12 12])
title('SWS')
subplot(3,2,4)
imagesc(b,axE,cer(cT==2,:)),colorbar
caxis([-2.5 6])
xlim([-12 12])
title('REM')
subplot(3,2,5)
errorbar(b,mean(cesw(cT==2,:)),sem(cesw(cT==2,:)),'k'),colorbar
caxis([-2.5 6])
xlim([-12 12])
ylabel('Average z-score')
xlabel('time-lag (ms')

title('SWS-Wake')

subplot(3,2,6)
errorbar(b,mean(cerw(cT==2,:)),sem(cerw(cT==2,:)),'k'),colorbar
caxis([-2.5 6])
xlim([-12 12])
title('Wake-REM')
ylabel('Average z-score')
xlabel('time-lag (ms')
title('REM-Wake')

if 0
    epswrite([figDir 'HDCell_MonoSyn_ToPYR.eps'],'size','screen')
end

axE = (1:sum(cT==1));

figure(3),clf
subplot(3,2,1)
imagesc(b,axE,ce(cT==1,:)),colorbar
caxis([-2.5 8])
xlim([-12 12])
ix9 = find(ixExcBo' &hdC & cIx(:,2)==9);
title('PostSyn INT, All states')
ylabel('Thal-Pos Pair #')

subplot(3,2,2)
imagesc(b,axE,cew(cT==1,:)),colorbar
xlim([-12 12])
caxis([-2.5 8])
title('Wake')
subplot(3,2,3)
imagesc(b,axE,ces(cT==1,:)),colorbar
xlim([-12 12])
caxis([-2.5 8])
title('SWS')
subplot(3,2,4)
imagesc(b,axE,cer(cT==1,:)),colorbar
xlim([-12 12])
caxis([-2.5 8])
title('REM')
subplot(3,2,5)
errorbar(b,mean(cesw(cT==1,:)),sem(cesw(cT==1,:)),'k'),cb=colorbar;set(cb,'visible','off')
caxis([-2.5 6])
xlim([-12 12])
ylim([-1.2 1.2])
ylabel('Average z-score')
xlabel('time-lag (ms')

title('SWS-Wake')
subplot(3,2,6)
errorbar(b,mean(cerw(cT==1,:)),sem(cerw(cT==1,:)),'k'),cb=colorbar;set(cb,'visible','off')
caxis([-2.5 6])
xlim([-12 12])
ylim([-1.2 1.2])
ylabel('Average z-score')
xlabel('time-lag (ms')
title('REM-Wake')



figure(4),clf
plot(b,mean(ce(cT==1,:)),'b')
hold on
plot(b,mean(ce(cT==2,:)),'r')
ylabel('Average z-score')
xlabel('time-lag (ms')
legend('INT','PYR')


postS = 19;

ixPrePair = find(ixExcBo' & cIx(:,2)==postS);
ixPreCell = cIx(ixPrePair,1);

ixPrePair = ixPrePair(find(hdCellStats(ixPreCell,end)));
ixPreCell = ixPreCell(find(hdCellStats(ixPreCell,end)));

str = mean(xcSAll(b>1.5&b<3,ixPrePair));
str = str/max(str);
figure(100),clf
    polar_ap(B,AngHisto{postS}/hdCellStats(postS,4),'r')
    hold on
    for ii=1:length(ixPreCell)
        polar(B,str(ii)*AngHisto{ixPreCell(ii)}/hdCellStats(ixPreCell(ii),4))
        %polar(B,AngHisto{ixPreCell(ii)}/hdCellStats(ixPreCell(ii),4))
        
    end









