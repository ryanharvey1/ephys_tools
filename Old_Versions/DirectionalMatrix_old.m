function [ Overall_DirectionalityIndex ] = DirectionalMatrix(tfile,event,timestamps,Overall_DirectionalityIndex,SmoothRateMap_Right,SmoothRateMap_Left)
%DirectionalMatrix Summary of this function goes here
%   Detailed explanation goes here
% NORMALIZATION
for k=1:size(SmoothRateMap_Right,1)
    for kk=1:size(SmoothRateMap_Right,2)
        SmoothRateMap_Right_Norm(k,kk) = (SmoothRateMap_Right(k,kk)-min(SmoothRateMap_Right(k,:)))/(range(SmoothRateMap_Right(k,:)));
        SmoothRateMap_Left_Norm(k,kk) = (SmoothRateMap_Left(k,kk)-min(SmoothRateMap_Left(k,:)))/(range(SmoothRateMap_Left(k,:)));
    end
end

% ARRANGE BINS
kkk=1;
for k=1:size(SmoothRateMap_Right_Norm,2)
    index=SmoothRateMap_Right_Norm(:,k)==1;
    for kk=1:size(SmoothRateMap_Right_Norm,1)
        if index(kk,1)==1
            SmoothRateMap_Right_arranged(kkk,:)=SmoothRateMap_Right_Norm(kk,:);
            SmoothRateMap_Left_arranged(kkk,:)=SmoothRateMap_Left_Norm(kk,:);
            kkk=kkk+1;
            continue
        end
    end
end

% Overall_DirectionalityIndex
Overall_DirectionalityIndex=(sum(Overall_DirectionalityIndex(:,1))-...
    sum(Overall_DirectionalityIndex(:,2)))/(sum(Overall_DirectionalityIndex(:,1))+sum(Overall_DirectionalityIndex(:,2)));

% PLOT LEFT VS RIGHT SMOOTHED RATEMAPS
scrsz=get(groot,'ScreenSize');
figure2=figure('OuterPosition',[1,scrsz(4)/2,scrsz(3)/2,scrsz(4)/2]); subplot(2,1,1),h = pcolor(SmoothRateMap_Right_arranged);
shading interp
colormap jet
axis off
hold on
colorbar
box off
set(h, 'EdgeColor', 'none');
set(gca,'YDir','reverse');
title('Smoothed Rate Map Right');

figure2; subplot(2,1,2),h = pcolor(SmoothRateMap_Left_arranged);
shading interp
colormap jet
axis off
hold on
colorbar
box off
set(h, 'EdgeColor', 'none');
set(gca,'YDir','reverse');
title('Smoothed Rate Map Left');

[filepath, ~] = fileparts(tfile);
if timestamps==1; % SAVE WITH S# IF EVENTS EXIST
    saveas(figure2,[filepath sprintf('S%d',event) '_SmoothedRatePlot.tiff']);
    save([filepath sprintf('S%d',event) '_Data.mat']);
else
    saveas(figure2,[filepath '_SmoothedRatePlot.tiff']);
    save([filepath '_Data.mat']);
end
close all
end



