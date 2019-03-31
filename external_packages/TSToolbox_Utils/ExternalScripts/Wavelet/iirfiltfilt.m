% IIRFILTFILT   filter a signal using a double N pole butterworth.
%
% call          Y = iirfiltfilt(X,Fs,flp,fhp,N,band)
%
%               XF = IIRFILTFILT(X, FS, FLP)
%               XF = IIRFILTFILT(X, FS, [], FHP)
%               XF = IIRFILTFILT(X, FS, FLP, FHP)
%               XF = IIRFILTFILT(...,N) determined number of poles (default - 3)
%               XF = IIRFILTFILT(...,BAND) {1} - bandpass; 0 generates bandstop
%
%               if X is a matrix, its columns are filtered.
%
% see also      FILTFILT, IIRFILT
%
% NOTES         1. handling of matrix data is not perfect
%               2. magnitude response is double the expected

% 13-jan-06 ES

% revisions
% 28-oct-11 made autonomous of the iirlp2hp/iirlp2bp/ routines

function y = iirfiltfilt(x,Fs,flp,fhp,N,band)

if nargin<6 || isempty(band), band = 1; end
if nargin<5 || isempty(N), N = 3; end
if any(size(x)==1), x = x(:); end
Y = [];

if nargin<3 || (isempty(flp) && isempty(fhp))
    error('at least one band limit mustbe specified')
elseif nargin==3 || isempty(fhp),  % lowpass
    wc = flp/Fs*2;
    [b a] = butter(N,wc,'low');
elseif nargin>=4 && isempty(flp),
    wc = fhp/Fs*2;
    [b a] = butter(N,wc,'high');
elseif nargin>=4 && band==1
    wc = [flp fhp]/Fs*2;
    [b a] = butter(N,wc );
elseif nargin>=4
    wc = [flp fhp]/Fs*2;
    [b a] = butter(N,wc,'stop' );
end

% prepare filter
len = size(x,1);
nb = length(b);
na = length(a);
nfilt = max(nb,na);
nfact = 3*(nfilt-1);
if (len<=nfact)
    error('Data must have length more than 3 times filter order.');
end
if nb < nfilt, b(nfilt)=0; end
if na < nfilt, a(nfilt)=0; end
rows = [1:nfilt-1  2:nfilt-1  1:nfilt-2];
cols = [ones(1,nfilt-1) 2:nfilt-1  2:nfilt-1];
data = [1+a(2) a(3:nfilt) ones(1,nfilt-2)  -ones(1,nfilt-2)];
sp = sparse(rows,cols,data);
zi = sp \ ( b(2:nfilt).' - a(2:nfilt).'*b(1) );
% flip edges (matrix wise)
pad = ones( nfact, 1 ); 
y = [2*pad*x(1,:)-x((nfact+1):-1:2,:);x;2*pad*x(len,:)-x((len-1):-1:len-nfact,:)];
% filter, reverse data, filter again, and reverse data again
y = filter(b,a,y,zi*y(1),1);
y = y(length(y):-1:1,:);
y = filter(b,a,y,zi*y(1),1);
y = y(length(y):-1:1,:);
% remove edges
y([1:nfact len+nfact+(1:nfact)],:) = [];

return
