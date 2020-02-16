function ncs2dat(varargin)
% ncs2dat: calls python file ncs2dat.py
% to run, cd to session folder and run this function

channels=1:length(dir('*.ncs'));

targetfile=pwd;

% create target file name (session folder name)
pathsplt=strsplit(targetfile,filesep);
targetfile=fullfile(targetfile,pathsplt{end});

% locate ncs2dat.py (should be within Utils)
py_loc=which('ncs2dat.py');

Args = struct('basename','CSC','zeropad',0,'pyFileLoc',['python ', py_loc],'fileExt','.ncs');
Args.flags = {'zeropad'};
[Args,~] = getOptArgs(varargin,Args,'remove',{});


files = [];
for ch = channels
    if Args.zeropad
        files  = cat(2,files,sprintf('%s%.02d%s ',Args.basename,ch,Args.fileExt));
        
    else
        files = cat(2,files,sprintf('%s%d%s ',Args.basename,ch,Args.fileExt));
    end
   
end

cmd = [Args.pyFileLoc ' ' files targetfile '.dat'];

[status,cmdout] = dos(cmd,'-echo');

if status == 0 
    display(['saved ' pwd filesep targetfile])
else
    display('ERROR ERROR ERROR ERROR')
%     save('ncs2datError.mat','cmdout')
    
end

function [ArgStruct,args]=getOptArgs(args,ArgStruct,varargin)
% getOptArgs Helper function for parsing varargin. 
%    [ArgStruct,ARGS] = getOptArgs(ARGS,ArgStruct,VARARGIN) parses 
%    ARGS for parameter-value pairs and enters them into ArgStruct,
%    which is a structure containing named arguments with default
%    values. Any numeric arguments at the beginning of ARGS that are
%    not matched to a parameter string are stored in the field 
%    NumericArguments of ArgStruct. The ARGS cell array, which may 
%    be modified via optional input arguments, is also returned.
%
%    The optional input arguments are:
%    'flags' - followed by a cell array indicating which arguments are 
%              flagtype arguments, i.e. arguments that don't require a 
%              value (the value will be set to 1 if it is present). If
%              this argument is not specified, but there is a field in
%              ArgStruct named 'flags', that will be used instead. If
%              this argument is present, it will be added to ArgStruct
%              if it is not already present, and will replace the flags
%              field if present.
%    'aliases' - followed by a cell array indicating aliases, which
%                can be used to map one argument-name to several 
%                argstruct fields.
%    'remove' - followed by a cell array indicating which arguments 
%               in ARGS are to be removed if present.
%    'subtract' - following by a cell array indicating which arguments
%                 in ARGS will have their values subtracted by 1, with
%                 a minimum value of 0.
%    'shortcuts' - followed by a cell array with the name of the 
%                  shortcut in the first column and a cell array in the
%                  second column with the full arguments.
%    'stopOnError' - if present, indicates that the function should
%                    stop execution when it finds unknown arguments.
%    example usage: 
%    --------------
%    function parseargtest(varargin)
%
%    %define the acceptable named arguments and assign default values
%    Args=struct('Holdaxis',0, ...
%           'SpacingVertical',0.05,'SpacingHorizontal',0.05, ...
%           'PaddingLeft',0,'PaddingRight',0,'PaddingTop',0, ...
%           'PaddingBottom',0,'MarginLeft',.1,'MarginRight',.1, ...
%           'MarginTop',.1,'MarginBottom',.1,'rows',[],'cols',[], ...
%           'RedoLevels',0,'SaveLevels',0,'Dir','','flags',{'Holdaxis'}); 
%
%    % The capital letters define abrreviations.  
%    % Eg. parseargtest('spacingvertical',0) is equivalent to  
%    % parseargtest('sv',0) 
%    % fill the arg-struct with values entered by the user
%    [Args,varargin]=getOptArgs(varargin,Args,'flags',{'Holdaxis'}, ... 
%        'aliases',{'Spacing' {'sh','sv'}; ...
%                   'Padding' {'pl','pr','pt','pb'}; ...
%                   'Margin' {'ml','mr','mt','mb'}}, ...
%        'shortcuts',{'redo',{'RedoLevels',1};'save',{'SaveLevels',1}}, ...
%        'subtract',{'RedoLevels','SaveLevels'},...
%        'remove',{'Dir'}, ...
%        'stopOnError');
%
%    disp(Args)

% Based on parseArgs.m from Matlab Central by Aslak Grinsted 2003

Aliases={};
FlagTypeParams='';
RemoveParams='';
SubtractParams='';
ShortcutParams='';
ShortCuts='';
% flag to indicate if we want to stop execution with an error or 
% continue
bError = 0;

% get length of varargin
nArgs = nargin - 2;
i = 1;

% check if ArgStruct has a field called flags
if(isfield(ArgStruct,'flags'))
    FlagTypeParams = lower(ArgStruct.flags);
end

% look for optional arguments
while(i <= nArgs)
    arg = varargin{i};	
    if ischar(arg)
        arg = lower(arg);
        switch arg
            case('flags')
                % make argument case insensitive
                FlagTypeParams = lower(varargin{i+1});
                i = i + 2;
                % set flags field in ArgStruct
                ArgStruct.flags = FlagTypeParams;
            case('aliases')
                Aliases = varargin{i+1};
                i = i + 2;
            case('remove')
                RemoveParams = lower(varargin{i+1});
                i = i + 2;
            case('subtract')
                SubtractParams = lower(varargin{i+1});
                i = i + 2;
            case('shortcuts')
                ShortCutParams = varargin{i+1};
                ShortCuts = lower(strvcat(ShortCutParams{:,1}));
                i = i + 2;
            case('stoponerror')
                bError = 1;
                i = i + 1;
            otherwise
                i = i + 1;
        end
    else
        i = i + 1;
    end
end

%
% Skip "numeric" arguments preceeding first param,value pair
%
% get number of arguments in args
numargs = size(args,2);
% start at 1 instead of 0 and subtract 1 later so we can actually 
% index into args
NumArgCount=1;
while (NumArgCount<=numargs)&(~ischar(args{NumArgCount}))
    NumArgCount=NumArgCount+1;
end
if(NumArgCount>1)
    ArgStruct.NumericArguments = {args{1:(NumArgCount-1)}};
else
    ArgStruct.NumericArguments = [];
end

%
% Make an accepted fieldname matrix (case insensitive)
%
Fnames=fieldnames(ArgStruct);
for i=1:length(Fnames)
    name=lower(Fnames{i,1});
    Fnames{i,2}=name; %col2=lower
    % find characters that are uppercase
    AbbrevIdx=find(Fnames{i,1}~=name);
    % store uppercase letters as abbreviations
    %col3=abreviation letters (those that are uppercase in the ArgStruct) 
    % e.g. SpacingHoriz->sh
    Fnames{i,3}=[name(AbbrevIdx) ' '];
    %the space prevents strvcat from removing empty lines
    %Does this parameter have a value? (e.g. not flagtype)
    Fnames{i,4}=isempty(strmatch(Fnames{i,2},FlagTypeParams,'exact'));
end
% convert to character matrices with 1 string in each row
FnamesFull=strvcat(Fnames{:,2});
FnamesAbbr=strvcat(Fnames{:,3});

nAliases = size(Aliases,1);
if nAliases>0  
    for i=1:nAliases
        name=lower(Aliases{i,2});
        % get length of aliases
        naliases = size(Aliases{i,2},2);
        % find indices corresponding to aliases
        FieldIdx = [];
        for j = 1:naliases
            idx = strmatch(name{j},FnamesAbbr,'exact');
            if isempty(idx)
                idx = strmatch(name{j},FnamesFull,'exact');
            end
            if ~isempty(idx)
                FieldIdx = [FieldIdx; idx];
            end
        end
        Aliases{i,2}=FieldIdx;
        name = lower(Aliases{i,1});
        % find characters that are uppercase
        AbbrevIdx=find(Aliases{i,1}~=name);
        % the space prevents strvcat from removing empty lines
        Aliases{i,3}=[name(AbbrevIdx) ' '];
        %dont need the name in uppercase anymore for aliases
        Aliases{i,1}=name;
    end
    %Append aliases to the end of FnamesFull and FnamesAbbr
    FnamesFull=strvcat(FnamesFull,strvcat(Aliases{:,1})); 
    FnamesAbbr=strvcat(FnamesAbbr,strvcat(Aliases{:,3}));
end

%--------------get parameters--------------------
l = NumArgCount;
al = length(args);
while (l<=al)
    a=args{l};
    if ischar(a) 
        if ~isempty(a)
            a=lower(a);
            %try abbreviations (must be exact)
            FieldIdx=strmatch(a,FnamesAbbr,'exact'); 
            if isempty(FieldIdx) 
                FieldIdx=strmatch(a,FnamesFull,'exact'); 
            end
            if FieldIdx>length(Fnames) %then it's an alias type.
                FieldIdx=Aliases{FieldIdx-length(Fnames),2}; 
            end
            
            if isempty(FieldIdx) 
                % check if it is a shortcut
                ShortCutIdx = strmatch(a,ShortCuts,'exact');
                if(~isempty(ShortCutIdx))
                    % replace shortcut with full string in args
                    if(l==al)
                        % if we are already at the end of the cell array
                        args = {args{1:(l-1)} ShortCutParams{ShortCutIdx,2}{:}};
                    elseif(l==1)
                        % if we are at the beginning of the cell array
                        args = {ShortCutParams{ShortCutIdx,2}{:} args{l+1:end}};
                    else
                        args = {args{1:(l-1)} ShortCutParams{ShortCutIdx,2}{:} args{l+1:end}};
                    end
                    al = length(args);
                elseif bError
                    error(['Unknown named parameter: ' a])
                else
                    l = l + 1;
                end
            else
                for curField=FieldIdx' %if it is an alias it could be more than one.
                    if (Fnames{curField,4})
                        val=args{l+1};
                    else
                        %parameter is of flag type and is set (1=true)....
                        val=1; 
                    end
                    ArgStruct.(Fnames{curField,1})=val;
                    if(strmatch(a,SubtractParams,'exact'))
                        % make sure value does not go below 0
                        args{l+1} = max([0 args{l+1} - 1]);
                    end
                end
                if(strmatch(a,RemoveParams,'exact'))
                    if(Fnames{FieldIdx(1),4})
                        % parameter with value so remove both argument and 
                        % value
                        [args,al] = removeargs(args,l,2);
                    else
                        % flag type so just remove argument
                        [args,al] = removeargs(args,l,1);
                    end
                else
                    % if flag type then go to next argument but if 
                    % param/value, then skip next argument
                    l=l+1+Fnames{FieldIdx(1),4}; 
                end
            end
        else 
            l=l+1;
        end
    else
        if bError
            error(['Expected a named parameter: ' num2str(a)])
        else
            l = l + 1;
        end
    end
end

