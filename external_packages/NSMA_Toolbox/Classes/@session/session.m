function ses = session(s)
% SESSION session class constructor:
%  ses = session()   or
%  ses = session(s)
%
%  Class holding all info associated with a recording session
%  designed to aid in analysis routines that loop over list of 
%  sessions and need to build t-file lists and corresponding S-matrices 
%  (=: cellarrays of spiketimes of each cell) per session epoch.
%  
%  s ... either a session class or a struct with fields (omitted fields are set by default to an empty string):
%        s.animal    ... animal (rat) ID  (string)
%        s.name      ... session name (string)
%        s.group     ... session group (string or cell array of strings)  (e.g treatment, cell type, recording system number,...);
%                        a session can be member of several session groups listed in the cell array of groups, e.g.
%                        {'CPP','pyramidal','SystemD'}; one can then build session lists based on membersip in the animal and
%                        group categories.
%        s.split     ... 0/1 (false/true) flag, describing if session is split into separate tfiles for each epoch;
%                        for split sessions, the epochs field holds tfileglob strings describing the tfiles for each epoch
%                        instead of the epoch start/end timestamps
%                        
%        s.tfiledir  ... path to tfile directory (string or cellarray of strings); in case of a cellarray of strings all
%                        direcories in the list are searched for tfiles sequentially.
%                        Subdirectories of the dirs in the list are NOT searched!
%        s.tfileglob ... glob string or cell array of stings defining all tfiles in tfiledir for this session 
%                        (string with a wildcard char *, e.g. 'rat8972_ses01_cell*.t);
%                        in case of a cellarray of strings, all strings are globbed sequentially and appended to the 
%                        tfile list.
%                        s.tfileglob is only used for non-split sessions! split session require the tfileglob strings
%                        to be given per epoch (see s.epochs). 
%        s.sameTTfmt ... format of the string identifying the tetrode number in the tfile name:
%                        e.g.:   '6401_Dec_8_tt%%' where the chars up to the %% are the same for all cells and %%
%                                 denotes the position of the tetrode number in the string.
%        s.epochs    ... an epochs struct, one field for each epoch, holding a vector with epoch start and end timestamps
%                        in case of continuous sessions, e.g.:
%                           s.epochs.sleep1 = [ts1,ts2]
%                           s.epochs.maze   = [ts3,ts4]
%                           s.epochs.sleep2 = [ts5,ts6]
%                        or a tfileglob string or a cell array of glob strings e.g :
%                           s.epochs.sleep1 = '6401_Dec_8_tt*_sleep1.*'
%                           s.epochs.maze = {'6401_Dec_8_tt*_maze1.*','6401_Dec_8_tt*_run1.*'}
%                           s.epochs.sleep2 = '6401_Dec_8_tt*_sleep2.*'
%                        The S-matrices for each epoch are built from the tfile lists resulting from these glob strings.
%                        The number of epochs and their names are user definable, but should be consitent within an analysis
%                        sequence.
%        s.other     ... a user defined struct holding extra info associated with this session that an analysis routine might want:
%                        e.g. filenames of eeg files, tracker position files, ripple times file, etc...:
%                        for example:
%                           s.other.posfile = 'VT1.pascci'
%                           s.other.posfiledir = 'c:\data\myana\rat6812\ses01'
%
%   The following fields are computed in constructor from the info given above:
%
%       s.tfilenames ... a cell array of strings with the names of the all the tfiles found in with
%                          s.tfileglob (in split=0 case) or s.epochs (in split=1 case) for all epochs.
%                          Only the filename.ext is given in the list, NOT the full path.
%       s.tfilefullnames ... same as above, however the full path is returned (e.g. C:\data\mytfile.t)
%
% Version 0.1.0
% PL Jan 2003

%defining the class fields and their defaults
ses.animal = '';
ses.name = '';
ses.group = '';
ses.split = '';
ses.tfiledir = '';
ses.tfileglob = '';
ses.sameTTfmt = '';
ses.epochs = '';
ses.other = '';
%
% computed fields
ses.tfilenames = '';
ses.tfilefullnames = '';

if nargin == 0
    % default constructor
    ses = class(ses,'session');
elseif isa(s,'session')
    % copy constructor
    ses = s;
elseif isa(s,'struct')
    % construct from input struct
    % populate all available fields given in input s
    fn = fieldnames(ses);
    for i=1:length(fn)
        if isfield(s,fn{i})
            val = getfield(s,fn{i});
            ses = setfield(ses,fn{i},val);
        end
    end
    
    % get the tfilenames:
    if exist(ses.tfiledir,'dir')
        pushdir(ses.tfiledir);
    else
        popdir
        error(['tfiledir: ', ses.tfiledir,' does not exist!']);    
    end
    if ses.split
        % accumulate tfiles for all epochs
        epochs = fieldnames(ses.epochs);
        fnlist = {};
        for i=1:length(epochs)
            tfileglob = getfield(ses.epochs,epochs{i});
            fnl = sort(findfiles(tfileglob,'CheckSubdirs',0));
            if isempty(fnl)
                popdir
                error([' there are no tfiles matching ' tfileglob ' in dir: ' ses.tfiledir]);
            end
            if i>1 & length(fnlist(:,1)) ~= length(fnl)
                popdir
                error(['tfiles in dir ' ses.tfiledir ...
                   ' do not match! You need to have the same number of tfiles for each epoch!']);
            end
            fnlist = [fnlist, fnl];
        end
        ses.tfilefullnames = fnlist;
        [nrows,ncols] = size(ses.tfilefullnames);
        ses.tfilenames = cell(nrows,ncols);
        for ic = 1:ncols
            ses.tfilenames(:,ic) = RemoveDirectoryFromFilename(ses.tfilefullnames(:,ic));
        end
    else
        ses.tfilefullnames = sort(FindFiles(ses.tfileglob,'CheckSubdirs',0));
        if ~isempty(ses.tfilefullnames) 
            ses.tfilenames = RemoveDirectoryFromFilename(ses.tfilefullnames);
        else
            ses.tfilenames = {};
        end
    end
    popdir;

    % make class from struct
    ses = class(ses,'session');

end