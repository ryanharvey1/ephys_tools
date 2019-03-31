function CalibrateLEDorientation_Moser()

[~,fbasename,~] = fileparts(pwd);

[X,Y,ang,wakeEp] = LoadPosition_Moser(fbasename);
ang0 = atan2(diff(Data(Y)),diff(Data(X)));
ang1 = Data(ang);
angOffset = CircularMean(ang1(2:end)-ang0);
if angOffset > pi
    angOffset = angOffset-2*pi;
end
figure(1),clf
plot(mod(ang0,2*pi),mod(ang1(1:end-1),2*pi),'k.')
title(num2str(angOffset))

figure(2),clf
plot(mod(ang0,2*pi),mod(ang1(1:end-1)-angOffset,2*pi),'k.')

info = {'Angle difference from Red front / Blue rear condition'};
SaveAnalysis(pwd,'GeneralInfo',{angOffset},{'angOffset'},info);
