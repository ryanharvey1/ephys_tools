ReadCR.cpp and the binary ReadCR.mexw32 are updated version (May 19,2008) with a bug fix found by Jean-Marc:
 The start_ts was computed incorrectly to return the record AFTER the one containing the requested start timesamp. This is now fixes.
The current old .mex32 dll is still in the IO dir and locked (accessed by somebody). The new version temporarily resides here till
I can move it to its proper place.

Your code using ReadCR.mex32 (called interally e.g. by ReadCR_tsd.m and others) should still work fine without changes. The differince is
that now one more record is returned at the beginning of your eeg-snippet (and that 
record now CONTAINS the sample point with the requested start timestamp)

PL May 19, 2008

(PS: if I forget to move this - plz someone remind me??? ThnX)