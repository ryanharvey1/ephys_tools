function ep = LoadEpoch(fbasename,epname)

%ep = LoadEpoch(fbasename,epname)

%Adrien Peyrache, 2012

Fs = 1250;

%try

    if exist([fbasename '-states.mat'],'file');
        load([fbasename '-states.mat'])
        states = states(:);
        t = [0:length(states)-1]';
        switch lower(epname)
            case 'sws'
                epIx = states==2 | states==3;
            case 'rem'
                epIx = states==5;
            case 'is'
                epIx = states==4;
            case 'wake'
                epIx = states==1;
        end
        epIx(1)=0;
        epIx(end)=0;
        ep = tsd(t,double(epIx));
        ep = thresholdIntervals(ep,0.5);
         switch lower(epname)
            case 'sws'
                ep = mergeCloseIntervals(ep,5);
                ep = dropShortIntervals(ep,60);
            case 'rem'
                ep = mergeCloseIntervals(ep,5);
                ep = dropShortIntervals(ep,10);
            case 'wake'
                ep = mergeCloseIntervals(ep,20);
                ep = dropShortIntervals(ep,60);
         end
                
    elseif exist([fbasename '.sts.' epname],'file')
        ep = load([fbasename '.sts.' epname]);
        ep = intervalSet(ep(:,1)/Fs,ep(:,2)/Fs);
    else
        warning(['LoadEpoch: no states or STS file to load for ' fbasename])
        ep = intervalSet([],[]);
    end

% catch
%     warning(lasterr);
%     keyboard
% end
