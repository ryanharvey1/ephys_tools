% Input/output functions (including batch handling) for FMAToolbox.
%
% NOTE
%
%   Although individual data files can be handled by these 'low-level' functions,
%   it is generally easier to use <a href="matlab:help SetCurrentSession">SetCurrentSession</a> and related <a href="matlab:help Data">Get...</a> functions
%   instead.
%
% Batch file handling
%
%   StartBatch           - Start a new batch job.
%   GetBatch             - Get batch job output.
%   BatchInfo            - Get batch job information.
%   CancelBatch          - Cancel batch job.
%   CleanBatches         - Delete completed batch jobs from memory.
%
% Creating new event files
%
%   NewEvents            - Create events structure.
%   SaveEvents           - Save events to evt file.
%   SaveRippleEvents     - Save ripple events to evt file.
%
% Low-level binary input/output
%
%   LoadBinary           - Load data from a binary file.
%   ResampleBinary       - Resample binary file (e.g. dat->lfp).
%   SaveBinary           - Save data to binary file.
%
% Low-level data input/output
%
%   LoadEvents           - Load events from an evt file.
%   LoadParameters       - Load parameters from an xml file.
%   LoadPositions        - Load position data from a pos file.
%   LoadSpikeTimes       - Load spike times from a res file.
%   LoadSpikeFeatures    - Load spike features from file.
%   LoadSpikeWaveforms   - Load spike waveforms from a spk file.
%