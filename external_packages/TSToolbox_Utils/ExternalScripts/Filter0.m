% y = Filter0(b, x)
%
% filters x with a fir filter so it has zero phase, i.e. shifts the
% filtered signal to the right half of the length of b.
%
% for now it zero pads the original signal
% later it might also do reflecton boundary conditions.
%
% be careful about the order of b!
% for even filters it is not exact  (change of Anton)
% - tired that even filterss dont' work


function y = Filter0(b, x)

if size(x,1) == 1
	x = x(:);
end

% if mod(length(b),2)~=1
% 	error('filter order should be odd');
% end
if mod(length(b),2)~=1
    shift = length(b)/2;
else
    shift = (length(b)-1)/2;
end

x=double(x);
[y0 z] = filter(b,1,x);

y = [y0(shift+1:end,:) ; z(1:shift,:)];


