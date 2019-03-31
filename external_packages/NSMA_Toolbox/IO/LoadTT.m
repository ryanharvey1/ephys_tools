function TT = LoadTT(fn)

% LoadTT  Reads an NSMA TT file
%
% TT = LoadTT(fn)
%
% INPUTS:
%   fn = .tt file
% OUTPUTS:
%   TT = tsd structure where data = nSpikes x nSamplesPerSpike x nTrodes
%
% Uses mex file LoadTT0(fn) to do the main read
%
% ADR 1998, version L5.1, last modified '03 by PL

% status PROMOTED
% PL 2003: added selection of TT files either by 'TT*' or '*.tt' or '*.ntt'


[fdir,fname,fext] = fileparts(fn);
if strcmpi(fname(1:2),'TT') | strcmpi(fext,'.tt') | strcmpi(fext,'.ntt')
   [t, wv] = ReadTT(fn);
elseif strcmpi(fname(1:2),'ST')
   [t, wv] = ReadST(fn);
elseif strcmpi(fname(1:2),'SE')
   [t, wv] = ReadSE(fn);
else
   warndlg('Only filenames starting with TT,ST or SE are recognized by LoadTT!!!','LoadTT Warning!');
   return;
end%if
%t(end) = [];
%wv(end,:,:)=[];
TT = tsd(floor(t),wv);












%--------------------------------------------------------------------------------
function oldLoadTT

nTrodes=4;		% 4 channels for tetrode
nSamplesPerSpike=32;	% number of samples per spike waveform per channel
maxSpikesToRead = inf;  % maximum spikes to read
debugON = 0;

Extract_varargin;	

%---- Open file
[fp, msg] = fopen(fn,'rb','b');
if (fp == -1);   error(msg); end
[tmpdn, tmpfn] = fileparts(fn);
TTWindowName = ['Load TT (' tmpfn ')'];

%---- Read Header
H = ReadHeader(fp);
HSkip = ftell(fp);

%--- Read Timestamps
if (debugON); disp('Reading timestamps...'); end
Ttt = fread(fp, maxSpikesToRead, 'uint32', 2*nSamplesPerSpike*nTrodes); 
nSpikes = length(Ttt);
if (debugON); disp('   Done.'); end

%--- Read V
if (debugON); disp('Reading waveforms...'); end
V = zeros(nSpikes, nTrodes, nSamplesPerSpike);
fseek(fp, HSkip, 'bof');
for iS = 1:nSpikes
   DisplayProgress(iS, nSpikes, 'Title', TTWindowName);
   fseek(fp, 4, 'cof');
   TMP = fread(fp, [nTrodes nSamplesPerSpike], 'int16');
   if (size(TMP) == [nTrodes nSamplesPerSpike])
      V(iS,:,:) = TMP;
   else
      warning('TT file ends early, ignoring last sample');
   end
end
if (debugON); disp('   Done.'); end

fclose(fp);

%--- generate T
TT = tsd(Ttt, V);

return;


