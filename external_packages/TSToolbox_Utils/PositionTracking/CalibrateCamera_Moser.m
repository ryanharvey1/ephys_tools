function CalibrateCamera()

Fs = 1250/32;

[dumy fbasename dumy] = fileparts(pwd);
[whl,t,GoodRanges] = LoadPosition(fbasename);
load('Analysis/BehavEpochs','wakeEp')

X = Restrict(tsd(t*10000,whl(:,1)),wakeEp);
Y = Restrict(tsd(t*10000,whl(:,2)),wakeEp);

figure(1),clf
plot(Data(X),Data(Y))

done=0;
while ~done
    fprintf('Draw a line on the figure\n')
    gg = ginput(2);
    lpixel = norm(gg(1,:)-gg(2,:));
    lcms = input('What is the actual size of this line (in cms)','s');
    bin2cms = str2num(lcms)/lpixel;
    fprintf('Ratio cms/pixel is: %f\n',bin2cms)
    answer = input('Satisfied? [Y/N]','s');
    if strcmp(answer,'Y')
        done=1;
    elseif strcmp(answer,'N')
        fprintf('OK, then let''s restart!\n')
    else
        fprintf('Don''t understand the answer, so we restart\n')
    end
end
 
info = {'Bins to cms ratio'};
SaveAnalysis(pwd,'GeneralInfo',{bin2cms},{'bin2cms'},info);
