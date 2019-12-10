function X = randvonmises (N, MU, K)
%  RANDVONMISES (N, MU, K)
%    Generates N random samples from a Von Mises distribution with mean MU
%    and concentration K.
%
%    Ref: Best, D. and Fisher, N. (1979). Applied Statistics, 24, 152-157.
%
%   --> www.paulbays.com

if K==0, X = (rand(N,1)*2-1)*pi; return; end

a = 1 + (1 + 4 * (K^2))^0.5;
b = (a - (2 * a)^0.5)/(2 * K);
r = (1 + b^2)/(2 * b);

obs = 1;

while (obs <= N)
    z = cos(pi * rand);
    f = (1 + r * z)/(r + z);
    c = K * (r - f);
    U = rand;
    if (c * (2 - c) - U > 0) | (log(c/U) + 1 - c >= 0)
        X(obs,1) = wrap(sign(rand - 0.5) * acos(f) + MU);
        obs = obs + 1;
    end
end