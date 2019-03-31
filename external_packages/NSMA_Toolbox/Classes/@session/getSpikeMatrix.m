function S = getSpikeMatrix(ses,epoch)
% S = getSpikeMatrix(ses,epoch)
% 
% for a given session object, load all available tfiles into a Spike Matrix (S-Matrix) restricted 
% to the given epoch.
% The order of the S-Matrix corresponds to the order of the t-files returned 
% by getTFileList(ses,epoch,fullpath).
% 
% INPUTS:
%   ses  ... the current session object
%   epoch .. string; must be a valid field name of the session.epochs struct.
%
% OUTPUT:
%   S   ... Spike Matrix (S-Matrix), a cell array of ts-objects holding the timestamps 
%           (in 0.1 msec units) of each tfile.
%
% PL Feb. 2003

if nargin == 1
    epoch = 'any';
end

if ses.split & ~isfield(ses.epochs,epoch)
    disp('For split session you must provide a valid epoch string as second input argument!');
    error(['Input argument: ' epoch ' is not a valid field in the session.epochs struct!']);
end

if exist(ses.tfiledir,'dir')
    pushdir(ses.tfiledir);
else
    popdir;
    error(['tfiledir: ', ses.tfiledir,' does not exist!']);    
end
if ses.split
    fnlist = getTFileList(ses,epoch,1);
    if isempty(fnlist)
        popdir;
        error([' there are no tfiles matching ' tfileglob ' in dir: ' ses.tfiledir]);
    end
    quiet = 1;
    S = LoadSpikes(fnlist,quiet);
    popdir;
else
    fnlist = getTFileList(ses,'any',1);
    quiet = 1;
    S = LoadSpikes(fnlist,quiet);
    ts_start_end = getfield(ses.epochs,epoch);
    for i = 1:length(S)
        S{i} = Restrict(S{i},ts_start_end(1),ts_start_end(2));
    end
    popdir;
end
