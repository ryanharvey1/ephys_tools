% AdaptiveBinning
% based on (Skaggs, 1996)

% ********TESTING**********

% Time per bin
occ =[1,2,5,7,6,5,40,60,20,3,4,10,11,4,1];
% Spikes per bin
Smatrix=[5,0,1,0,0,0,10,20,19,0,0,1,1,1,1];
% scaling factor



% MAIN LOOP

binned=[]; Smatrixbin=Smatrix(1); occbin=occ(1); binNum=1; Scaling=100;
for i=1:length(occ)
    if Smatrixbin>Scaling/(occbin^2*(binNum/2)^2)
        FR=Smatrixbin/occbin; binned=[binned,FR];
        if i+1>length(occ); break; end
        binNum=1; Smatrixbin=Smatrix(i+1); occbin=occ(i+1);
        continue
    end
    if i+1>length(occ); 
        binNum=2; Smatrixbin=Smatrixbin+Smatrix(i-1); occbin=occbin+occ(i-1);
        for ii=1:length(occ)
            if Smatrixbin>Scaling/(occbin^2*(binNum/2)^2)
                FR=Smatrixbin/occbin; binned=[binned,FR]; 
                break
            end
            binNum=binNum+1; Smatrixbin=Smatrixbin+Smatrix(i-1); occbin=occbin+occ(i-1);
        end
        break 
    end
    binNum=binNum+1; Smatrixbin=Smatrixbin+Smatrix(i+1); occbin=occbin+occ(i+1);
end

% RESIZE
binned=imresize(binned,[1,length(Smatrix)],'nearest');
% COMPARE

figure(1); subplot(4,1,1); pcolor([(Smatrix./occ);(Smatrix./occ)])
 colormap jet
title('Normal Binning')
subplot(4,1,2); pcolor([(Smatrix./occ);(Smatrix./occ)])
 colormap jet; shading interp
title('Normal Binning')

subplot(4,1,3); pcolor([binned;binned])
colormap jet
title('Adaptive Binning')
subplot(4,1,4); pcolor([binned;binned])
colormap jet; shading interp
title('Adaptive Binning')


