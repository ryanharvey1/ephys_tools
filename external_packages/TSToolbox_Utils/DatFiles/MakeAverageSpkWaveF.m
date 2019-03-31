function MakeAverageSpkWaveF(fbasename)

% This function has to be launched from root data folder and will load the
% different dat files in each of the fbasename-xx folders

rootname = [fbasename filesep fbasename];
xmlFname = [rootname '.xml'];

analysisDir = [fbasename filesep 'Analysis'];

if ~exist(analysisDir,'dir')
   mkdir(analysisDir);   
end

if ~exist(xmlFname,'file')
    error(['XML file ''' xmlfname ''' does not exist'])
end
warning off
xmldata = xml_load([rootname '.xml']);
warning on
segnames = xmldata.multiFileProcessing.files;
xmldata = LoadXml([rootname '.xml']);
nChannels = xmldata.nChannels;
nShanks = length(xmldata.SpkGrps);

waveM_AllShanks = cell(nShanks,1);
waveStd_AllShanks = cell(nShanks,1);

segfname = segnames(1).fileBaseName;
segfname = [segfname filesep segfname];
try
    for sh=1:nShanks
        fprintf('Shank %d...\n',sh);
%        if exist([rootname '.clu.' num2str(sh)],'file') && exist([rootname '.res.' num2str(sh)],'file')
%            clu = load([rootname '.clu.' num2str(sh)]);
%            res = load([rootname '.res.' num2str(sh)]);

        if exist([segfname '.clu.' num2str(sh)],'file') && exist([segfname '.res.' num2str(sh)],'file')
            clu = load([segfname '.clu.' num2str(sh)]);
            res = load([segfname '.res.' num2str(sh)]);

            clu = clu(2:end);
            cluIx = unique(clu);
            cluIx = cluIx(cluIx>1);
            cluIx = cluIx(:)';
            
            chanIx = xmldata.SpkGrps(sh).Channels;
            chanIx = chanIx(:)';
            
            waveM = [];
            waveStd = [];
            
            for ch=1:length(chanIx)
                fprintf('  Channel %d\n',chanIx(ch));
                dat = [];
                for seg=1:1 %length(segnames); Let's do it on the first one...
                    fname = segnames(seg).fileBaseName;
                    fprintf('\t%s...\n',fname);
                    fname = [fname filesep fname];
                    dat = [dat;LoadBinary([fname '.dat'],'channels',chanIx(ch)+1,'nChannels',nChannels,'frequency',20000)];
                end
                
                for c=1:length(cluIx)
                    fprintf('\tCell %d\n',c);
                    t = res(clu==cluIx(c));
                    [EegSegAv, EegSegStd] = TriggeredAvM_Data(dat,t,[1 2],20000);
                    if isempty(waveM)
                        waveM = zeros(length(cluIx),length(chanIx),length(EegSegAv));
                        waveStd = zeros(length(cluIx),length(chanIx),length(EegSegAv));
                    end
                    waveM(c,ch,:) = EegSegAv;
                    waveStd(c,ch,:) = EegSegStd;
                end                      
                
            end

        end
        waveM_AllShanks{sh} =  waveM;
        waveStd_AllShanks{sh} =  waveStd;
    end

    save([analysisDir filesep 'SpkWaveForms.mat'],'waveM_AllShanks','waveStd_AllShanks');
    
catch
    warning(lasterr);
    keyboard
end
