% [CleanWhl GoodRanges] = CleanWhl(Whl, StretchLen, JumpSize);
%
% "cleans up" a wheel file by interpolating missing stretches
% up to StretchLen long (default 20), for which the endpoints 
% don't differ by more than JumpSize, (default 30).
%
% also returns the ranges where the whl file is valid (in .whl units)
% GoodRanges which gives start and end samples of the good ranges
% (so Whl(GoodRanges) has no -1 values).
%
% if there are any very high derivatives left over, it warns you

function [cWhl, GoodRanges] = CleanWhl(Whl, StretchLen, JumpSize)

if nargin<2
	StretchLen = 20;
end

if nargin<3
    JumpSize = 30;
end

nWhl = size(Whl,1);

% interpolate missing values or large jumps
BigJump = abs(diff(Whl(:,1)))>10 | abs(diff(Whl(:,2)))>10;
Good = find(Whl(:,1)>-1 & ~([BigJump;0] | [0;BigJump]));

if length(Good)<2
    cWhl = -ones(size(Whl));
else
    cWhl = round(interp1(Good, Whl(Good,:), 1:nWhl, 'linear', -1));
end

% find missing stretches
dGood = [-(Whl(1,1)==-1) ; diff(Whl(:,1)>-1)];
BadStart = find(dGood<0);
BadEnd = find(dGood>0)-1;
if BadEnd(1)<BadStart(1)
    BadEnd(1)=[];
end
% if last point is bad, need to finish off BadEnd
if Whl(end,1)==-1
	BadEnd = [BadEnd; nWhl];
end

% find ranges to chuck
% jump size ...
if length(BadStart>0)
    StartInd = clip(BadStart-1, 1, nWhl); % StartInd and EndInd give the 
    EndInd = clip(BadEnd+1, 1, nWhl);     % points you are interpolating between

	ToChuck = find(BadEnd-BadStart>=StretchLen ...
		| abs(Whl(StartInd,1)-Whl(EndInd,1)) > JumpSize ...
		| abs(Whl(StartInd,2)-Whl(EndInd,2)) > JumpSize);
	
	% chuck em
	for i=ToChuck(:)'
       	cWhl(BadStart(i):BadEnd(i),:) = -1;
	end
end

if 0 % OLD VERSION (BUG?)
% % now find good ranges
% dcGood = [-(Whl(1,1)==1) ; diff(cWhl(:,1)>-1)];
% GoodStart = find(dcGood>0);
% GoodEnd = find(dcGood<0)-1;
% % if last point is good, need to finish GoodEnd
% if cWhl(end,1)>-1
%     GoodEnd = [GoodEnd; nWhl];
% end
% GoodRanges = [GoodStart, GoodEnd];
else
    dcGood = diff([0; cWhl(:,1)>-1; 0]);
    GoodStart = find(dcGood>0);
    GoodEnd = find(dcGood<0)-1;
    GoodRanges = [GoodStart, GoodEnd];
end


% delete singletons
if length(GoodStart>0)
    Singletons = find(GoodStart==GoodEnd);
    cWhl(GoodStart(Singletons),:) = -1;
    GoodRanges(Singletons,:) = [];
end

%keyboard


return
% find points that are not -1 and no big jumps (the rat can't move 3cm in 25 ms)
BigJump = abs(diff(Whl(:,1)))>10 | abs(diff(Whl(:,2)))>10;
Good = find(Whl(:,1)>-1 & ~([BigJump;0] | [0;BigJump]));


% interpolate missing values
cWhl = round(interp1(Good, Whl(Good,:), 1:size(Whl,1), 'linear', -1));

% now overwrite bad bits

dGood = [-(Whl(1,1)==-1) ; diff(Whl(:,1)>-1)];

BadStart = find(dGood<0);
BadEnd = find(dGood>0)-1;

% if last point is bad, need to finish off BadEnd
if Whl(end,1)==-1
	BadEnd = [BadEnd; size(Whl,1)];
end

TooLong = find(BadEnd-BadStart>StretchLen);
for i=TooLong(:)'
	cWhl(BadStart(i):BadEnd(i),:) = -1;
end

WhlGood = [0 ; Whl(:,1)~=-1; 0];
dWhlGood = diff(WhlGood);
GoodStart = find(dWhlGood==1);
GoodEnd = find(dWhlGood==-1)-1;
LongEnough = (GoodEnd-GoodStart)>StretchLen;
GoodRanges = [GoodStart(LongEnough), GoodEnd(LongEnough)]; % in .whl units
