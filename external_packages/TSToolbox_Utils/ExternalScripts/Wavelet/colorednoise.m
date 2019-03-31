% colorednoise      pink/red/blue/violet noise
%
% y = colorednoise( N, color )
%
% N - number of samples to be returned in row vector
% color - {'pink'}, i.e. 1/f. alternatives include 1/f^2 (red), f (blue), and f^2 (violet).
% y - row vector of colored noise samples
%
% The function generates a sequence of colored noise samples.
% Pink noise has equal energy in all octaves (or similar log bundles) of frequency.
% In terms of power at a constant bandwidth:
%   pink noise falls off at 3 dB per octave.
%   red noise falls off at 6 dB per octave.
%   blue noise increase in at 3 dB per octave.
%   violet noise increase in at 6 dB per octave.

% 30-sep-13 ES
% based on scripts by Hristo Zhivomirov

function y = colorednoise( N , color )

nargs = nargin;
if nargs < 1 || isempty( N )
    N = 1000; 
end
if nargs < 2 || isempty( color )
    color = 'white'; 
end

% define the length of the vector
% ensure that the M is even
if rem(N,2)
    M = N+1;
else
    M = N;
end

% generate white noise with sigma = 1, mu = 0
x = randn( M, 1 );

% prepare a vector for multiplication
switch color
    case 'pink' % 1/f
        NumUniquePts = M/2 + 1;
        n = 1:NumUniquePts;
        n = sqrt(n);
    case 'red'  % 1/(f^2)
        NumUniquePts = M/2 + 1;
        n = 1:NumUniquePts;
    case 'blue' % f multiplication
        NumUniquePts = M/2 + 1;
        n = 1:NumUniquePts;
        n = sqrt(n);
    case 'violet' % f^2 multiplication
        NumUniquePts = M/2 + 1;
        n = 1:NumUniquePts;
    otherwise % white
        y = x;
        return
end

% FFT
X = fft(x);

% multiplicate the left half of the spectrum so the power spectral density
switch color
    case { 'pink', 'red' }
        % pink:
        % is inversely proportional to the frequency by factor 1/f, i.e. the
        % amplitudes are inversely proportional to 1/sqrt(f)

        % red:
        % is inversely proportional to the frequency by factor 1/(f^2), i.e. the
        % amplitudes are inversely proportional to 1/f
        X(1:NumUniquePts) = X(1:NumUniquePts)./n( : );
        
        
    case { 'blue', 'violet' }
        % blue:
        % is proportional to the frequency by factor f, i.e. the
        % amplitudes are proportional to sqrt(f)
        
        % violet:
        % is proportional to the frequency by factor f^2, i.e. the
        % amplitudes are proportional to f
        X(1:NumUniquePts) = X(1:NumUniquePts).*n( : );
end

% prepare a right half of the spectrum - a copy of the left one,
% except the DC component and Nyquist frequency - they are unique
X(NumUniquePts+1:M) = real(X(M/2:-1:2)) -1i*imag(X(M/2:-1:2));

% IFFT
y = ifft(X);

% prepare output vector y
y = real(y(1:N));

% normalize
y = y./max(abs(y));

return

