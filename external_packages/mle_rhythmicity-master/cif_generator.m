function [ cif, cif_int ] = cif_generator( funname )
%CIF_GENERATOR Generate function handles for rhythm conditional intensity 
% This helper function generates function handles for fitting rhythmic
% distributions of lags. Using functional programming in this way greatly
% speeds up evaluation during runtime.
%
% INPUT
% funname: String name of the CIF handles to be generated. Can be:
%   flat - Exponential decay with no oscillation
%   pure - Pure sinusoidal decay
%   noskip - Decaying rhythm with no skipping
%   full - Decaying rhythm with skipping
%   noskip_rise - Noskip with exponential rise
%   full_rise - Decaying rhythm with skipping and exponential rise
%
% OUTPUT
%   cif - Function handle for the Conditional Intensity Function (CIF).
%       with the form @(t,...). t is the lags at which
%       to evaluate the CIF, and it is followed by the appropriate
%       parameters
%   cif_int - The function handle for the definite integral of the CIF,
%       with the form @(t,...). 
%
% Copyright 2015-2016 Trustees of Boston University
% All rights reserved.
%
% This file is part of mle_rhythmicity revision 2.0. The last committed
% version of the previous revision is the SHA starting with 93862ac...
%
% This code has been freely distributed by the authors under the BSD
% license (http://opensource.org/licenses/BSD2-Clause). If used or
% modified, we would appreciate if you cited our papers:
%
% Climer JR, DiTullio R, Newman EL, Hasselmo ME, Eden UT. (2014),
% Examination of rhythmicity of extracellularly recorded neurons in the
% entorhinal cortex. Hippocampus, 25:460-473. doi: 10.1002/hipo.22383.
%
% Hinman et al., Multiple Running Speed Signals in Medial Entorhinal
% Cortex, Neuron (2016). http://dx.doi.org/10.1016/j.neuron.2016.06.027

switch funname
    case 'flat'% Non-rhythmic
        cif = @(t,tau,b)(1-b)*exp(-t*10^-tau)+b;
        cif_int = @(t,tau,b)10^tau*(1-b)*(1-exp(-t*10^-tau))+b*t;
    case 'flat_rise'% Not-rhythmic with an exponential rise
        cif = @(t,tau,b,v)((1-b)*exp(-t*10^-tau)+b).*(1-2.^(-t/v));
        cif_int = @(t,tau,b,v)b*t+(b-1)*10^tau*exp(-t*10^-tau)+b*2.^(-t/v)*v/log(2)+(1-b)*2^(-t/v)*exp(-t*10^-tau)*v/(10^-tau*v+log(2));
        cif_int = @(t,tau,b,v)cif_int(t,tau,b,v)-cif_int(0,tau,b,v);
    case 'pure'% Pure sinusoid
        cif = @(t,f)cos(2*pi*f*t)+1;
        cif_int = @(t,f)sin(2*pi*f*t)/(2*pi*f)+t;
    case 'noskip'% Non-skipping
        cif = @(t,tau,b,c,f,r)(1-b).*exp(-t.*10.^-tau).*(r.*exp(-t.*10.^-c).*...
           cos(2*pi*f.*t)...
           +1)+b;
       cif_int = @(t,tau,b,c,f,r)b*t+(1-b).*10.^tau.*(1-exp(-t*10.^-tau))+r.*(1-b).*...
           exp(-t*(10.^-c+10.^-tau)).*((10.^-c+10.^-tau).*(exp(t.*(10.^-c+10.^-tau))-cos(2*pi*f*t))+2*pi*f.*sin(2*pi*f*t))./...
           (4*pi^2*f.^2+(10.^-c+10.^-tau).^2);
    case 'full'
        fixb = @(b)max(b,realmin);% b>0
        cif = @(t,tau,b,c,f,s,r)(1-b).*exp(-t.*10.^-tau).*(r.*exp(-t.*10.^-c).*...
           ((2+2*sqrt(1-s)-s).*cos(2*pi*f.*t)+4*s.*cos(pi*f.*t)+2-2*sqrt(1-s)-3*s)/4 ...
           +1)+b;
       cif = @(t,tau,b,c,f,s,r)cif(t,tau,fixb(b),c,f,s,r);
       cif_int = @skipping_integral;
    case 'noskip_rise'
       cif = cif_generator('noskip');
       cif = @(t,tau,b,c,f,v,r)cif(t,tau,b,c,f,r).*(1-2.^(-t./v));
       
       cif_int = @noskip_rise_integral;
    case 'full_rise'
        cif = cif_generator('full');
        cif = @(t,tau,b,c,f,s,v,r)cif(t,tau,b,c,f,s,r).*(1-2.^(-t./v));
       cif_int = @full_skipping_integral;     
    otherwise
end

end

% Integral of the full distribution with skipping
function D = skipping_integral(t,tau,b,c,f,s,r)
    A = 10.^-c+10.^-tau;
    B = exp(t*A);
    b = max(b,realmin);
    
    D=b*t+...
        ...
        (1-b).*...
        (10.^tau).*...
        (1-exp(-t*10.^-tau))+...
        (r.*(1-b)/4).*...
    ((2+2*sqrt(1-s)-s).*((B.^-1).*(A.*(B-cos(2*pi*t*f))+2*pi*f.*sin(2*pi*t*f))./...
    (4*pi^2*f.^2+A.^2))+...
    4*s.*((B.^-1).*(A.*(B-cos(pi*t*f))+pi*f.*sin(pi*t*f))./...
    (pi^2*f.^2+A.^2))+...
    (2-2*sqrt(1-s)-3*s).*(10.^(c+tau)).*(1-(B.^-1))./(10.^c+10.^tau));
end

function D = noskip_rise_integral(t,tau,b,c,f,v,r)
    A = 10.^-tau+10.^-c;
    B = A+log(2)./v;
    
    D = @(t,tau,b,c,f,v,r)b*t...
        +b*2.^(-t./v).*v/log(2)...
        -(1-b).*10.^tau.*exp(-t*10.^-tau)...
        -(b-1).*exp(-t*(log(2)./v+10.^-tau))/(log(2)./v+10.^-tau)...
        +(1-b).*r.*exp(-t*A).*(-A.*cos(2*pi*t*f)+2*pi*f.*sin(2*pi*t*f))./(A.^2+4*pi^2*f.^2)...
        +(b-1).*r.*exp(-t*B).*(-B.*cos(2*pi*t*f)+2*pi*f.*sin(2*pi*t*f))./(B.^2+4*pi^2*f.^2);
    D = D(t,tau,b,c,f,v,r)-D(0,tau,b,c,f,v,r);
end

function D = full_skipping_integral(t,tau,b,c,f,s,v,r)
    A = 10.^-tau+10.^-c;
    B = A+log(2)./v;

D = @(t,tau,b,c,f,s,v,r)b*t...
    +1/log(2)*b.*2.^(-t./v).*v...
    +(b-1).*exp(-t*10.^-tau).*10.^tau...
    +(1-b).*exp(-t.*(log(2)./v+10.^-tau))./(log(2)./v+10.^-tau)...
    +r.*(1-b).*(2*(sqrt(1-s)-1)+3*s)/4.*exp(-t*A)./A...
    +r.*(b-1).*(2*(sqrt(1-s)-1)+3*s)/4.*exp(-t*B)./B...
    +r.*s.*(1-b).*exp(-t*A).*(-A*cos(pi*t*f)+pi*f.*sin(pi*t*f))./(A.^2+pi^2*f.^2)...
    +r.*s.*(b-1).*exp(-t*B).*(-B*cos(pi*t*f)+pi*f.*sin(pi*t*f))./(B.^2+pi^2*f.^2)...
    + 0.25*(1-b).*(2*(1+sqrt(1-s))-s).*r.*exp(-t*A).*(-A.*cos(2*pi*t*f)+2*pi*f.*sin(2*pi*t*f))./(A.^2+4*pi^2*f.^2)...
    + 0.25*(b-1).*(2*(1+sqrt(1-s))-s).*r.*exp(-t*B).*(-B.*cos(2*pi*t*f)+2*pi*f.*sin(2*pi*t*f))./(B.^2+4*pi^2*f.^2);

D = D(t,tau,b,c,f,s,v,r)-D(0,tau,b,c,f,s,v,r);
end