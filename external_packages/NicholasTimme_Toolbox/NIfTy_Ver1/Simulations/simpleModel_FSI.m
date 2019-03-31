%%FS - Layer 5 rat visual cortex fast-spiking (FS) interneuron (fig8.27 from 2007 book)
function [mv,spkCt,spkLoc]=simpleModel_FSI(vInput,T,tau,PLOT,nzInput)

%% FS interneuron parameters 
CI=20; vrI=-55; vtI=-40; kI=3.5;
aI=0.2; bI=-20; cI=-45; dI=-55; %% Note, Chris originally had aI = 0.1, fig 8.27 uses aI = 0.2
vpeak=25;

%% Variables/Stim times for running as a script
% PLOT=1;
% pulseSize=550;
% T=100;tau=0.001;    % dt (comment for script)

n=round(T/tau);       % number of samples 
mv=vrI*ones(1,n); u=0*mv;

% vInput=[zeros(1,round(n*0.25)) ones(1,n*0.5).*pulseSize zeros(1,n*0.25)]; % Comment out if function
% nzInput=zeros(size(vInput));

%% Setting variables
spkCt=0;
spkLoc=[];
m=1;

for i=1:n-1                         % forward Euler method
    mv(i+1)=(mv(i)+tau*(kI*(mv(i)-vrI)*(mv(i)-vtI)-u(i)+vInput(i))/CI)+nzInput(i);
    % For FS neurons, include nonlinear U(v): U(v) = 0 when v<vb ; U(v) = 0.025(v-vb) when v>=vb (d=vb=-55)
    if (mv(i+1) < dI)
        u(i+1) = u(i) + tau*aI*(0-u(i));
    else
        u(i+1) = u(i) + tau*aI*((0.025*(mv(i)-dI).^3)-u(i));
    end
    %  Check if spike occurred and need to reset
    if mv(i+1)>=vpeak
        mv(i)=vpeak;
        mv(i+1)=cI;
        spkCt=spkCt+1;     % count number of spikes and location
        spkLoc=[spkLoc m];
    end
    m=m+1;  
end

if isempty(spkCt); spkCt=0; end;

if PLOT==1;
   %% Plotting tools
   figure;
   subplot(2,1,1);
   plotyy(tau*(1:n),mv,tau*(1:n),vInput);
   xlabel('time (ms)'); 
   ylabel('membrane potential (mV)');
   subplot(2,1,2);
   plot(mv,u);
   xlabel('membrane potential (mV)');
   ylabel('recovery variable, u');
end;
  
