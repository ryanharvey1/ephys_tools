%LoadSpikeWaveF - Load subset of spikes
%
%  USAGE
%
%    spk = LoadSpikeWaveF(filename,nChannels,nSamples,CluIx)
%
%    filename       file to read
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'filename'    
%     'nChannels'   
%     'nSamples'    
%     'CluIx'       boolean vector of the spike indices to be loaded
%    =========================================================================

% Adrien 2015

function data = LoadSpikeWaveF(filename,nChannels,nSamples,CluIx)

% Default values
start = 0;
precision = 'int16';
skip = 0;
duration = Inf;
frequency = 20000;
channels = [];

if ~exist(filename),
	error(['File ''' filename ''' not found.']);
end
f = fopen(filename,'r');
if f == -1,
	error(['Cannot read ' filename ' (insufficient access rights?).']);
end

nSamplesPerSpk = nChannels*nSamples;

data = [];
chunkSz = 10000;
nbChunks = floor(CluIx(end)/chunkSz);
data = zeros(nSamplesPerSpk*length(CluIx),1,'int16');
dIx=1;

for ii=1:nbChunks
    ix = CluIx>(ii-1)*chunkSz & CluIx<=ii*chunkSz;
    if length(ix)
    d = fread(f,nSamplesPerSpk*chunkSz,precision);
    d = reshape(d,nChannels,nSamples,[]);
    d = d(:,:,CluIx(ix)-chunkSz*(ii-1));
    d = d(:);
    data(dIx:dIx-1+sum(ix)*nSamplesPerSpk) = d;
    dIx = dIx+sum(ix)*nSamplesPerSpk;
    end
end

if any(CluIx>chunkSz*nbChunks)
    ix = CluIx>chunkSz*nbChunks;
    d = fread(f,nSamplesPerSpk*(CluIx(end)-chunkSz*nbChunks),precision);
    d = reshape(d,nChannels,nSamples,[]);
    d = d(:,:,CluIx(ix)-chunkSz*nbChunks);
    d = d(:);
    data(dIx:dIx-1+sum(ix)*nSamplesPerSpk) = d;
end

fclose(f);
data = reshape(data,nChannels,nSamples,[]);

