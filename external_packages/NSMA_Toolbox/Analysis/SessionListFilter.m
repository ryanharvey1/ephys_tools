function filteredSessList = SessionListFilter(sessList,varargin)
% filteredSessList = FilterSessionList(sessList,varargin)
%
% filters a session List according to arguments given in a sequence of pairs:
%   you can filter according to the 'animal', 'name', or 'group' fields of the session object.
%
%   vargin ... pairs of strings:  'field '  , 'value'
%                          e.g.:  'animal'  , '6408'
%                                 'group'   , 'CPP'  
%                                  etc.
%   The value part of a pair can be a cell array of value strings. In this case the session
%   is included if ANY one string of the values matches any one of the corresponding session field.
%   This is equivalent to a logical OR.
%   
%   If multiple pairs are given, they are treated as a logical AND.  The same filename can
%   appear several times in the list.
%              
%   PL Feb. 2003

fSL = {};
for isess = 1:length(sessList)
    ignore = 0;
    sess= sessList{isess};
    for iarg = 1:2:length(varargin)
        field = varargin{iarg};
        value = varargin{iarg+1};
        if isempty(strmatch(field,{'animal','name','group'}))
            error([field ' is not an allowed field for SessionListFilter']);
        end
        if ~iscell(value)
            value = {value};   % make a one element cell array
        end
        add = 0;
        for iv = 1:length(value)
            if ~isempty(strmatch(upper(value{iv}),upper(getfield(struct(sess),field))))
                add = add+1;   % logical OR between values in same cell array
            end
        end
        if ~add
           ignore = 1;    % logical AND between pairs of field,value
        end
    end
    
    if ~ignore
       fSL{end+1} = sess;    
    end
end

filteredSessList = fSL;