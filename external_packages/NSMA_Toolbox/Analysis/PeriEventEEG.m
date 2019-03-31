function  [PeriEEG,EEGts,avgPeriEEG,stdPeriEEG] = PeriEventEEG(eegtsd, eventts, halfwindowsize)
% 
% [PeriEEG,EEGts,avgPeriEEG,stdPeriEEG] = PeriEventEEG(eegtsd, eventts, halfwindowsize)
%
% Peri-Event-Time Histogram (PETH) of an EEG tsd  with respect to a list of trigger-event timestamps T.
% The window size around T (at 0) can be asymmetric and ranges from [T - tlow_msec, T + thigh_msec].
% The number of bins depends on the sample frequency in the EEG file and the length of your window.
%
%  INPUT:
%     eegtsd ..... tsd of filtered EEG
%     eventts .... ts object with timestamps of events
%     halfwindowsize .. Delta t before and after the event timestamp (msec)
%
% OUTPUT:
%   PeriEEG  .......  nEvt*nSamp matrix with EEG samples centered at each of Nevt event timstamp
%   EEGts    ....     1 * nSamp vector with timestamps of bin centers (in msec units!)
%   avgPeriEEG .....  1 * nSamp matrix with columnwise average of above
%   stdPeriEEG .....  1 * nSamp matrix with columnwise ERROR OF THE MEAN of above
%
%  PL 2001
%  updated PL Aug. 2009

evts = Range(eventts);

% check that all events are within eeg limits
ix1 = evts < StartTime(eegtsd);
if any(ix1)
    warning('%d events before StartTime of eeg are ignored!',sum(ix1));
end
ix2 = evts > EndTime(eegtsd);
if any(ix2)
    warning('%d events after EndTime of eeg are ignored!',sum(ix2));
end
ix12 = ix1 | ix2;
evts(ix12) = [];  % remove out-of-bounds events
    

dt = halfwindowsize * 10;  % convert windowsize to timestamp units
lowlimit = evts - dt;
upperlimit = evts + dt;
% create a list of eeg snippets centered around events
nevts = length(evts);
eeglist = cell(1,nevts);
for i = 1:nevts
   eeglist{i} = Restrict(eegtsd,lowlimit(i),upperlimit(i));  % list of eeg snippets
end % for

% find max length and bounds of all snippets
ts_lo = +Inf;
ts_hi = -Inf;
maxlength = 0;  
for i = 1:nevts
   tsList = Range(eeglist{i},'ts');
   if isempty(tsList)
       error('Event %d is out of range of eeg timestamps. Make sure ALL input events are within EEG range!',i)
   end
   if tsList(1)-evts(i) < ts_lo
       ts_lo = tsList(1)-evts(i);
   end
   if tsList(end)-evts(i) > ts_hi
       ts_hi = tsList(end)-evts(i);
   end
   ll = length(tsList);
   if ll > maxlength
      maxlength = ll;
   end %
end%for

% align and collect all snippets into an array
dt = (ts_hi-ts_lo)/maxlength;       % binsize of PETH
EEGts = (ts_lo:dt:ts_hi);           % get bin low edges
EEGts = EEGts(1:end-1) + dt/2;      % disard last bin edge and shift to bin centers  
[tsmin,ix00] = min(abs(EEGts));     % find index of 0 bin in PeriEEG matrix
PeriEEG = zeros(nevts, maxlength);  
for i = 1:nevts
   eegvect = Data(eeglist{i});
   dt = Range(eeglist{i})-evts(i);
   [tsmin, ix0] = min(abs(dt));     % find index of 0 bin in current eeg snippet
   % find position of current snippet wrt PeriEEG 0-bin
   nts = length(dt); 
   dix = ix00(1)-ix0(1);            % offset between bincenters 
   ix1 = max(1, 1-dix );            % find lower edge of aligned eegvect snippet
   ix2 = min(nts, maxlength-dix);   % find upper edge of aligned eegvect snippet
   ix11 = max(1,1+dix);             % lower edge in PeriEEG array
   ix22 = min(maxlength,nts+dix);   % upper edge in PeriEEG array
   if (ix2-ix1)~=(ix22-ix11)
       error('Peter screwed up');   % check that index boundaries give consistent overlap
   end
   PeriEEG(i,ix11:ix22) = eegvect(ix1:ix2)';   % copy aligned snippet to array
end%for

EEGts = EEGts/10;     % convert to msec  
avgPeriEEG = mean(PeriEEG);
stdPeriEEG = std(PeriEEG)/sqrt(length(avgPeriEEG));

   

