function tfList = getTFileList(ses,epoch,fullpath)
% tfList = getTFileList(ses,epoch,fullpath)
% 
% for a given session object, returns all available tfiles for a given epoch 
% in sorted cell array of strings. 
%
% INPUTS: 
%  ses  ....... the current session object
%  epoch ...... string, must be a valid field name of the session.epochs struct.
%               This argument can be omitted if you don't have split sessions.
%               If epoch is 'any' the first epoch name in ses.epochs is used. 
%  fullpath ... optional flag, set to 1 if you want the full pathnames prepended to the filenames,
%               omit or set to 0 otherwise.
%
% PL Feb. 2003

% check input
if nargin <= 2
    fullpath = 0;
end
if nargin == 1
    epoch = 'any';
end

if strcmpi(epoch,'any')
    % get first epoch name in ses.epochs
    ep = fieldnames(ses.epochs);
    epoch = ep{1};
end
if ses.split & ~isfield(ses.epochs,epoch)
    disp('For split session you must provide a valid epoch string as second input argument!');
    error(['Input argument: ' epoch ' is not a valid field in the session.epochs struct!']);
end

% cd to tfiledir
if exist(ses.tfiledir,'dir')
    pushdir(ses.tfiledir);
else
    popdir
    error(['tfiledir: ', ses.tfiledir,' does not exist!']);    
end

% retrieve filenames
if ses.split
    tfileglob = getfield(ses.epochs,epoch);
    fnlist = sort(findfiles(tfileglob,'CheckSubdirs',0));
    if isempty(fnlist)
        popdir;
        error([' there are no tfiles matching ' tfileglob ' in dir: ' ses.tfiledir]);
    end
    if fullpath
        tfList = fnlist;
    else
        tfList = RemoveDirectoryFromFilename(fnlist)';
    end
    popdir;
else
    if fullpath
        tfList = ses.tfilefullnames';
    else
        tfList = RemoveDirectoryFromFilename(ses.tfilefullnames)';
    end
    popdir;
end
