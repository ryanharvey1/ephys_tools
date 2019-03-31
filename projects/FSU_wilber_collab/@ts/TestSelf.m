function Errors = TestSelf

% ts/TestSelf  Runs tests of all ts methods on a fake ts object
%
% Errors = TestSelf
%
% INPUTS:
%       (none)
% OUTPUTS:
%       Errors - cell array of error messages resulting from failed tests

Errors = {};

t = 1:10;
x = ts(t);

% Data
if Data(x) ~= t
   Errors = [Errors, {'Data(x) failed.'}];
end

% EndTime
if EndTime(x) ~= 10
   Errors = [Errors, {'EndTime failed.'}];
end

% StartTime
if StartTime(x) ~= 1
   Errors = [Errors, {'StartTime failed.'}];
end

% Restrict
y = Restrict(x, 5, 6);
if Data(y) ~=  [5 6]
   Errors = [Errors {'Restrict(x, t0, t1) failed.'}];
end