function b = subsref(ses,index)
%SUBSREF Define field name indexing for session objects
%
% PL Feb. 2003

switch index.type
case '.'
    switch lower(index.subs)
    case 'animal'
        b = ses.animal;
    case 'name'
        b = ses.name;
    case 'group'
        b = ses.group;
    case 'split'
        b = ses.split;
    case 'tfiledir'
        b = ses.tfiledir;
    case 'tfileglob'
        b = ses.tfileglob;
    case 'samettfmt'
        b = ses.sameTTfmt;
    case 'epochs'
        b = ses.epochs;
    case 'tfilenames'
        b = ses.tfilenames;
    case 'tfilefullnames'
        b = ses.tfilefullnames;
    case 'other'
        b = ses.other;
    otherwise
        error('Invalid field name');
    end
case '{}'
    error('Cell array indexing not supported by session objects')
case '()'
    error('Array indexing not supported by session objects')
end