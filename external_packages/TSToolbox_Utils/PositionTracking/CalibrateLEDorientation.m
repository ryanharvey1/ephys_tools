function CalibrateLEDorientation(varargin)

% USAGE
% CalibrateLEDorientation(angOffset)
% 
% calibrate LED position relative to the expected position blue front / red Rear (to check...)
% inputs:
%   angOffset: offset value (optional)
  
  
%Parameters
checkOldOffset = 0;

load('Analysis/BehavEpochs.mat','wakeEp')

if checkOldOffset
    load('Analysis/GeneralInfo.mat','angOffset')
    fprintf('Old AngOffset: %f\n',angOffset);  
end

if ~isempty(varargin)
    angOffset = varargin{1};
    if ~isa(angOffset,'numeric')
        error('AngOffSet must be a numeric value')
    end
  
else
[~,fbasename,~] = extractfbasename(pwd);

[X,Y,ep,wstruct] = LoadPosition_Wrapper(fbasename);
[ang,ep] = HeadDirection_Wrapper(fbasename,0);
ep = intersect(ep,wakeEp);

ang = Restrict(ang,ep);
%X = Restrict(X,ep);
%Y = Restrict(Y,ep);
X = Restrict(tsd(wstruct.t,wstruct.whl(:,3)),ep);
Y = Restrict(tsd(wstruct.t,wstruct.whl(:,4)),ep);

ang0 = atan2(diff(Data(Y)),diff(Data(X)));
ang1 = Data(ang);
angOffset = CircularMean(ang1(2:end)-ang0);
if angOffset > pi
    angOffset = angOffset-2*pi;
end
figure(1),clf
plot(mod(ang0,2*pi),mod(ang1(1:end-1),2*pi),'k.')
title(num2str(angOffset))

 fprintf('New AngOffset: %f\n',angOffset); 
 
pause

end
info = {'Angle difference from Red front / Blue rear condition'};
SaveAnalysis(pwd,'GeneralInfo',{angOffset},{'angOffset'},info);
