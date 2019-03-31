function [rightEp,leftEp] = LoadLinearTrackRuns(fbasename)

[dir f1 f2] =  fileparts(fbasename);
fbasename = [f1 f2];
fbasename = [dir filesep fbasename filesep fbasename];
Fs = 1250/32;

try

    load([fbasename '.whlL.mat']);
    ix = ismember(Lwhl.cont,'right medium');
    t = Lwhl.Time{ix};
    rightEp = intervalSet(10000*t(:,1)/Fs,10000*t(:,2)/Fs);
    ix = ismember(Lwhl.cont,'left medium');
    t = Lwhl.Time{ix};
    leftEp = intervalSet(10000*t(:,1)/Fs,10000*t(:,2)/Fs);
    
catch
    warning(lasterr);
    keyboard
end
