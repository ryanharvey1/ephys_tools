function RemoveSpikes(ts, windowSize)

% RemoveSpikes(ts, windowSize)
%
% INPUTS
%     ts -- ts object containing timestamps to remove
%     windowSize -- spikes will be removed if a bad spike occurs within +/- windowsize of it
%
% Steps through all tt files in current directory and removes those
% spikes and writes a new tt file.
%
% ADR 1998
% version L4.0
% status UNDER CONSTRUCTION

StatusWarning('UNDER CONSTRUCTION', 'RemoveSpikes');

BadSpikes = Data(ts);

% First find all TT files 
% and then sort into order from largest to smallest
% for efficiency

ttfns = FindFiles('TT*.tt');
fsize = zeros(size(ttfns));
for iF = 1:length(ttfns)
  D = dir(ttfns{iF});
  fsize(iF) = D.bytes;
end
[dummy,order] = sort(fsize);
ttfns = ttfns(order(end:-1:1));

% For each TT file from largest to smallest 
% remove timestamps and rewrite it
for iF = 1:length(ttfns)
   [dn, fn, ext] = fileparts(ttfns{iF});
   newfn = [fn '_fixed' ext];  
   fprintf(2, '%s\n\t', fn);
   
   TT = LoadTT(ttfns{iF});
   T = Range(TT, 'ts');
   WV = Data(TT);
   nSpikes = length(T);
   
   for iS = 1:nSpikes
      DisplayProgress(iS, nSpikes, 'Title', fn);
      ix = binsearch(BadSpikes, T(iS));
      if abs(T(iS) - BadSpikes(ix)) < windowSize
         T(iS) = -1;
      end
      if (ix < length(BadSpikes) & abs(T(iS) - BadSpikes(i+1)) < windowSize
         T(iS) = -1;
      end       
   end
   
   f = find(T == -1);
   T(f) = [];
   WV(f, :, :) = [];   
   TT = tsd(T, WV);
   
   WriteTT(fullfile(dn, newfn), TT);
end
