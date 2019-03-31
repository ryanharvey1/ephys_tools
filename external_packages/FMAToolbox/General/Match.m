%Match - Replace values in one list with closest values in a second list.
%
%  USAGE
%
%    [matched,only1,only2] = Match(list1,list2,<options>)
%
%    list1          first list
%    list2          second list
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'match'       optional matching criterion (see below)
%     'error'       matched values for which the difference exceeds this
%                   threshold are set to NaN (default = inf)
%    =========================================================================
%
%  OUTPUT
%
%    matched        values in list2 that matched those of list1
%    only1          elements in list1 that were not matched by any element in list2
%    only2          elements in list2 that did not match any element in list1
%
%    only1 and only2 are lists of logical indices.
%
%  NOTE
%
%    Pairs are matched in the following way. Let list1 = [v(1)...v(n)] and
%    list2 = [w(1)...w(m)]. We are trying to match each element in list1 with
%    one of the elements of list2. Consider one element v(k) of list1:
%
%       w(1) <= ... <= w(j) <= v(k) <= w(j+1) <= ... <= w(m)
%
%    This will be matched with either w(j) or w(j+1), depending on the match
%    criterion:
%
%       1) match = 'down' (default) => w(j)
%       2) match = 'up' => w(j+1)
%
%    Note that this function is generally asymmetrical.

% Copyright (C) 2004-2010 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

function [matched,only1,only2] = Match(list1,list2,varargin)

% Default values
match = 'down';
err = inf;

% Check parameters
if nargin < 2 | mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help <a href="matlab:help Match">Match</a>'' for details).');
end
if ~isdvector(list1) | ~isdvector(list2),
	error('Incorrect sizes: ''list1'' and ''list2'' must be vectors (type ''help <a href="matlab:help Match">Match</a>'' for details).');
end
list1 = list1(:);
list2 = list2(:);

% Parse parameter list
for i = 1:2:length(varargin),
  if ~isa(varargin{i},'char'),
    error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help Match">Match</a>'' for details).']);
  end
  switch(lower(varargin{i})),

    case 'match',
		match = lower(varargin{i+1});
		if ~isstring(trim,'up','down'),
			error('Incorrect value for property ''match'' (type ''help <a href="matlab:help Match">Match</a>'' for details).');
      end

	case 'error',
		err = varargin{i+1};
		if ~isdscalar(err,'>=0'),
			error('Incorrect value for property ''error'' (type ''help <a href="matlab:help Match">Match</a>'' for details).');
		end

    otherwise,
      error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help Match">Match</a>'' for details).']);

  end
end

% Now compare [t1 t1 ... t1;t2 t2 ... t2;...;tn tn ... tn] with [l1 l2 ... lm;l1 l2 ... lm;...;l1 l2 ... lm]
where = repmat(list2,1,length(list1)) <= repmat(list1',length(list2),1);

% This yielded something like [1 1 1;1 1 1;...;0 1 1;...;0 0 1;...;0 0 0;...;0 0 0]
% where the transition from 1 to 0 in each column indicates the index of the last element
% in 'list1' inferior to the corresponding element in 'list2'. Now, find these transitions.
% (we append ones before and zeros after 'where' to handle the special cases where the first
% element in 'list1' inferior to the element of 'list2' is the last element of 'list1'
% or does not exist).
where = -diff([ones(1,size(where,2));where;zeros(1,size(where,2))]);

% Convert the resulting logical matrix into a list of subscript indices
n = size(where,1);
where = mod(find(where),n);
where(where==0) = n;
if strcmp(match,'down'),
	where = where-1;
end

% Special cases: when matching down (resp. up), values in list1 lesser (resp. greater) than all values in list2
% cannot be matched; set them to NaN
good = where>0 & where<n;
matched = nan(size(where));

% List corresponding values in list2
matched(good) = list2(where(good));

% Check error threshold
bad = abs(matched-list1) > err+eps; % eps is to compensate for inevitable numerical rounding errors such as 1.1-1.0>0.1
matched(bad) = nan;

% Output
only1 = isnan(matched);
matched = matched(~only1);
only2 = logical(ones(size(list2)));
only2(where) = 0;