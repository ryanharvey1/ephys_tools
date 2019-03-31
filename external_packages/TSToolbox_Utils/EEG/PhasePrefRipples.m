function A = PhasePref(A)

%Parameters
b = fir1(96,[0.1 0.3]);
%  cellId = 44;

A = getResource(A,'CellNames');
A = getResource(A,'SpikeData');
A = getResource(A,'HcTrace');

A = getResource(A,'Sleep1Epoch');
sleep1Epoch = sleep1Epoch{1};
A = getResource(A,'Sleep2Epoch');
sleep2Epoch = sleep2Epoch{1};
A = getResource(A,'MidRipS1SWS');
A = getResource(A,'MidRipS2SWS');

midRipS1SWS = Range(midRipS1SWS{hcTrace-4});
midRipS2SWS = Range(midRipS2SWS{hcTrace-4});
ripInt1 = intervalSet(midRipS1SWS - 250,midRipS1SWS + 250);
ripInt2 = intervalSet(midRipS2SWS - 250,midRipS2SWS + 250);

A = registerResource(A, 'Sleep1RipplesPhase', 'tsdArray', {[], []}, ...
    'sleep1RipplesPhase', ...
    ['spike ripples phase '],'mfile');


A = registerResource(A, 'Sleep2RipplesPhase', 'tsdArray', {[], []}, ...
    'sleep2RipplesPhase', ...
    ['spike ripples phase '],'mfile');


resdir = [parent_dir(A), filesep 'PhasePlot/Ripples'];
[p,ds,e] = fileparts(current_dir(A));

eegfname = [current_dir(A) filesep ds 'eeg' num2str(hcTrace) '.mat'];

if exist([eegfname '.gz'])
    display(['unzipping file ' eegfname]);
    eval(['!gunzip ' eegfname '.gz']);
end
load(eegfname)
%  display(['zipping file ' eegfname]);
%  eval(['!gzip ' eegfname]);

eval(['eegHc1 = Restrict(EEG' num2str(hcTrace) ',sleep1Epoch);']);
eval(['eegHc2 = Restrict(EEG' num2str(hcTrace) ',sleep2Epoch);']);
eval(['clear EEG' num2str(hcTrace) ';']);

dEegHc1 = Data(eegHc1);
dEegHc1 = filtfilt(b,1,dEegHc1);

dEegHc2 = Data(eegHc2);
dEegHc2 = filtfilt(b,1,dEegHc2);

eegHc1 = tsd(Range(eegHc1),dEegHc1);
eegHc2 = tsd(Range(eegHc2),dEegHc2);

phHc1 = ThetaPhase(S, eegHc1, Start(ripInt1), End(ripInt1));
phHc2 = ThetaPhase(S, eegHc2, Start(ripInt2), End(ripInt2));

%  
%  keyboard
%  PhaseMap(XS,YS,phHc{1});


for i = 1:length(S)
	
	if length(Data(phHc2{i}))>0

		fh1 = figure(1),clf
		rose(2*pi*Data(phHc1{i}));
		title('Phase Preference with pfc Theta')
		fh2 = figure(2),clf
		rose(2*pi*Data(phHc2{i}));
		title('Phase Preference with hc Theta')
		keyboard		

%  		saveas(fh1,[resdir filesep ds '_' cellnames{i} 'PhasePrefRipplesSleep1'],'png');
%  		saveas(fh2,[resdir filesep ds '_' cellnames{i} 'PhasePrefRipplesSleep2'],'png');

	end

end



%  pfcThetaPhase =  tsdArray(phPfc);
%  hcThetaPhase = tsdArray(phHc);

sleep1RipplesPhase = phHc1;
sleep2RipplesPhase = phHc2;

%  keyboard

A = saveAllResources(A);
