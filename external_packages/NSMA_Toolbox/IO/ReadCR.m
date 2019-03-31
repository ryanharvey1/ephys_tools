% Read a CR or CSC file from Cheetah (NT and unix Version 2.x)
% (MEX FILE) 
%
% [cr_ts, cr, sFreq] = ReadCR(fname, start_ts, end_ts);
%
% INPUTS: 
%   fname = full filename of Cheetah_NT CSC*.dat file 
%   start_ts, end_ts (optional) = start and end timestamps (0.1 msec units) of portion to read
%             (Note: the actual returned intervals may be shorter than the
%             requested intervals since the records are recorded in blocks
%             of 512 sample points. To be safe, add the ts-equivalent of
%             512 sample points on either end of the interval, depending on
%             your sampling rate!)
% OUTPUTS:
%   cr_ts = n x 1: timestamp of each CR Record in file
%   cr = n x 512   one CR record of 512 ADC samples
%   sFreq = 1 x 1 common Sampling Freqency of CR data (Hz) 
%
% IMPORTANT NOTE regarding CHEETAH SAMPLING FREQUENCIES:
%   older Cheetah versions report misleading sampling frequencies in the
%   EEG data files.  To be safe, ALWAYS compute the actual sampling rate from the data by matlab code 
%   such as:
%   eeg = ReadCR_tsd('CSC01.ncs')   % return an eeg tsd with interpolated
%                                   %    timestamps for each sample point
%   ts = Range(eeg,'ts');  % ts = array of timestamps in 0.1 msec units
%   sFreq = 10000/median(diff(ts));
%
% PL 1999, version 2.0, last modified 2/02 by SC
%
% cowen Oct 2000 - modified for partial load
% cowen Feb 2002 - modified to load old v1.? UNIX cheetah data.