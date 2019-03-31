% makegaussfir          for 1D smoothing
%
% gwin = makegaussfir( sd, Fs )
% 
% sd        Gaussian SD [sec]
% Fs        sampling rate [Hz]

% 26-dec-12 ES

function gwin = makegaussfir( sd, Fs, support )

n = 2 * floor( sd * Fs / 2 ) + 1; % make odd
if ~exist( 'support', 'var' )
    support = 6;
end
switch length( sd )
    case 1
        win = gausskernel( n, 0, support * n( 1 ) + 1, 1 ); % make w/ support of 3 SDs
    case 2
        win = gausskernel( n( 1 ), n( 2 ), support * n( 1 ) + 1, support * n( 2 ) + 1 ); % make w/ support of 3 SDs
end
gwin = win / sum( win( : ) );

return

% function K = gausskernel(sigmaX,sigmaY,N,M)
% x = -(N-1)/2:(N-1)/2;
% y = -(M-1)/2:(M-1)/2;
% if sigmaY == 0 && sigmaX == 0
%     K = zeros(M,N);
%     return;
% elseif sigmaY == 0
%     X = -inf*ones(M,N);
%     Y = X;
%     X(ceil(M/2),:) = x;
%     Y(ceil(M/2),:) = 0;
%     sigmaY = 1;
% elseif sigmaX == 0
%     Y = -inf*ones(M,N);
%     X = Y;
%     X(:,ceil(N/2)) = 0;
%     Y(:,ceil(N/2)) = y(:);
%     sigmaX = 1;
% else
%     X = repmat(x,M,1);
%     Y = repmat(y,N,1)';
% end
% K = 1/(2*pi*sigmaX*sigmaY)*exp(-(X.^2/2/sigmaX^2)-(Y.^2/2/sigmaY^2));
% return

% EOF

% for uncorrelated noise, boxcar filtering reduces the variance by sqrt( n )
% and gaussian filtering by sqrt( n ) * 1.96 (actually slightly less reduction d.t. edge effects)
nw = 5; r = randn( 1e6, 1 ); 
w = ones( nw, 1 ) / nw; rf = firfilt( r, w ); 
g = makegaussfir( nw, 1, 18 ); rg = firfilt( r, g ); 
[ std( r ), std( r ) / sqrt( nw ) std( rf ) std( rg ) std( rg ) * 1.96  ]
