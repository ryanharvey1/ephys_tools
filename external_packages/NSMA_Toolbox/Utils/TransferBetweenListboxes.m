function TransferBetweenListboxes

% TransferBetweenListboxes  Transfers values between two listboxes
%
% TransferBetweenListboxes
%
% INPUTS:
%       (none)
% OUTPUTS:
%       (none)
%
% Given two listboxes, each of which contains the other's handle
% as its UserData, and given this as each one's callback.  This
% then transfers values between the two by clicking on the values.
%
% ADR 1998, version L4.0, last modified '98 by ADR

% status PROMOTED


fromLB = gcbo;
toLB = get(fromLB, 'UserData');
valueFrom = get(fromLB, 'Value');

toString = get(toLB, 'String');
fromString = get(fromLB, 'String');

if isempty(toString)
  toString = fromString(valueFrom);
else
  toString = toString;
  toString(end+1) = fromString(valueFrom);
  toString = sortcell(toString);
end
fromString(valueFrom) = [];

set(fromLB, 'String', fromString);
set(toLB, 'String', toString);

set(fromLB, 'Value', max(1, min(length(fromString), valueFrom)));
