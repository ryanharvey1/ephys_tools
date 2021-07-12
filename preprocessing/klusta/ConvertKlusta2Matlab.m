function [fet, clu, spktimes]= ConvertKlusta2Matlab(shank,basepath,varargin)
% USAGE
%
% [fet clu spktimes wav]= ConvertKlusta2Matlab(shank,basepath,basename,wvformExtract,saveFiles,numCores)
%
% INPUTS
%
% shank - the shank number (as a number, not text) equalling the name of
%         the folder under basepath with the data of interst.  Default = 1.
% basepath - directory path to the main recording folder with .dat and .xml
%            as well as shank folders made by makeProbeMapKlusta2.m (default is
%            current directory matlab is pointed to)
% basename - shared file name of .dat and .xml (default is last part of
%            current directory path, ie most immediate folder name)
% saveFiles - binary (0 or 1) to choose whether to save clu, res, fet, wav(.spk)
%             files or just return them as outputs
%
% OUTPUTS
%
% fet - features of waveforms
% clu - cluster ID's for waveforms
% spktimes(res) - spike times of waveforms(in sample #'s, not seconds)
% wav- waveform shapes taken from raw .dat file
%
%Converts .kwik/kwx files from Klusta into klusters-compatible
% fet,res,clu,spk files.  Works on a single shank of a recording, assumes a
% 16bit .dat and an .xml file is present in "basepath" (home folder) and
% that they are named basename.dat and basename.xml.  Also assumes that
% subdirectories in that basepath are made for each shank with names
% specified by numbers (ie 1,2,3,4..8).  In each shank folder should be
% .kwik and .kwx files made by klusta with names as follows:
% basename_sh[shankumber].kwik/kwx.
%
%
% Brendon Watson 2016
% edited by David Tingley 1/2017
% edited by Aza and Antonio 7/2017
% edited by Laura Berkowitz for ephys_tools 2/2021

% Unpack 
p = inputParser;
addParameter(p,'saveFiles',1,@isnumeric)
p.parse(varargin{:})

saveFiles = p.Results.saveFiles;
[~,basename] = fileparts(basepath);

% Make file names
if saveFiles 
    cluname = fullfile(basepath,[basename '.clu.' num2str(shank)]);
    resname = fullfile(basepath,[basename '.res.' num2str(shank)]);
    fetname = fullfile(basepath,[basename '.fet.' num2str(shank)]);
end

% Start grabbing data
tkwik = fullfile(basepath,['klusta_',num2str(shank)],[basename '_sh' num2str(shank) '.kwik']);
tkwx = fullfile(basepath,['klusta_',num2str(shank)],[basename '_sh' num2str(shank) '.kwx']);
clu = h5read(tkwik,['/channel_groups/' num2str(shank) '/spikes/clusters/main']);

%% Getting spiketimes
spktimes = h5read(tkwik,['/channel_groups/' num2str(shank) '/spikes/time_samples']);

%% Spike features
fets = h5read(tkwx,['/channel_groups/' num2str(shank) '/features_masks']);
fets = double(squeeze(fets(1,:,:)));
fet = fets';
fetMultiplier = double(intmax('int32'))/max(abs(fets(:))); % masked klustakwik has small floats that need to be expanded before rounding below
fet = fet .* fetMultiplier;

%% writing to clu, res, fet
if saveFiles
        
    fid=fopen(cluname,'w');
    fprintf(fid,'%.0f\n',clu);
    fclose(fid);
    clear fid
    
    fid=fopen(resname,'w');
    fprintf(fid,'%.0f\n',spktimes);
    fclose(fid);
    clear fid
    %%
    SaveFetIn(fetname,fet);

    disp(['Shank ' num2str(shank) ' done'])
end

function SaveFetIn(FileName, Fet, BufSize)
if nargin<3 | isempty(BufSize)
    BufSize = inf;
end

nFeatures = size(Fet, 2);
formatstring = '%d';
for ii=2:nFeatures
    formatstring = [formatstring,'\t%d'];
end
formatstring = [formatstring,'\n'];

outputfile = fopen(FileName,'w');
fprintf(outputfile, '%d\n', nFeatures);

if isinf(BufSize)
    fprintf(outputfile,formatstring,round(Fet'));
else
    nBuf = floor(size(Fet,1)/BufSize)+1;
    
    for i=1:nBuf
        BufInd = [(i-1)*nBuf+1:min(i*nBuf,size(Fet,1))];
        fprintf(outputfile,formatstring,round(Fet(BufInd,:)'));
    end
end

