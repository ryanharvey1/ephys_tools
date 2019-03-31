function pos = FindDirOnPath(full_dirname)
%
%  checks if directory full_dirname is on the matlab path and returns its first 
%  position pos = 0,1,2,3,... etc  in the path list. Returns 0 if the directory is NOT
%  found on the path.
%
%  This is useful to check if a given directory is the top level dir on the matlab path and
%  allows a program to take according action

pp = path;        % pp is one long string of fulldirs separated by ';' (unix-style)

% make cell array of dirs from pp
plist = {};
tok = '0';          % start with any tok string
while ~isempty(tok)
    [tok,pp] = strtok(pp,';');
    if ~isempty(tok)
        plist{end+1} = lower(tok);
    end
end

% now find the reqested dir
indx = strmatch(lower(full_dirname),plist,'exact');   %indx is a possibly empty vector of list indices
if isempty(indx)
    pos = 0;
else
    pos = indx(1);    % pos of first occurence
end