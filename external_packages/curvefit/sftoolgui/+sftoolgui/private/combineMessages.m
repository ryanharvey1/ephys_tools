function newStr = combineMessages(strA, strB)
%COMBINEMESSAGES   Combine two messages.
%
%   NEWSTR = COMBINEMESSAGES(STRA, STRB) creates a new string where STRA and
%   STRB are combined and separated by a new line. If either is empty, the other
%   is returned with no new line.

%   Copyright 2009 The MathWorks, Inc.
%   $Revision: 1.1.10.1 $    $Date: 2009/03/09 19:07:25 $  

if isempty( strA )
    newStr = strB;
elseif isempty( strB )
    newStr = strA;
else
    newStr = sprintf( '%s\n%s', strA, strB );
end
