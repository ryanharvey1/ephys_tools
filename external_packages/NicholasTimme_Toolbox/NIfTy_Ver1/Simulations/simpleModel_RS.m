 %% Izhikevich simple neuron Chapter 8 Dynamical Systems in neuroscience
function [mv,spkCt,spkLoc]=simpleModel_RS(vInput,T,tau,PLOT,nzInput)
  
%%RS pyramidal neuron parameters
C=100; vr=-60; vt=-40; k=0.7;       
a=0.03; b=-2; c=-50; d=100;
vpeak=35;

%% Variables/Stim times for running as a script 
% PLOT=1;
% pulseSize=100;
% T=100;tau=0.001;    %dt (comment for script)

n=round(T/tau);     % number of samples 
mv=vr*ones(1,n); u=0*mv;

%vInput=[zeros(1,round(n*0.25)) ones(1,n*0.5).*pulseSize zeros(1,n*0.25)]; % Comment out if function

%% Setting variables
spkCt=0;
spkLoc=[];
m=1;

for i=1:n-1;
    mv(i+1)=(mv(i)+tau*(k*(mv(i)-vr)*(mv(i)-vt)-u(i)+vInput(i))/C)+nzInput(i);
    u(i+1)=u(i)+tau*a*(b*(mv(i)-vr)-u(i));
    if mv(i+1)>=vpeak;                % spike is fired!
        mv(i)=vpeak;                  % padding the spike amplitude
        mv(i+1)=c;                    % membrane voltage rest
        u(i+1)=u(i+1)+d;              % recovery variable update
        spkCt=spkCt+1;                % count number of spikes and location
        spkLoc=[spkLoc m];
    end;
    m=m+1;
end;

 
if PLOT==1;
    %% Plotting tools
    subplot(2,1,1);
    plotyy(tau*(1:n),mv,tau*(1:n),vInput);
    xlabel('time (ms)');
    ylabel('membrane potential (mV)');
    subplot(2,1,2);
    plot(mv,u);
    xlabel('membrane potential (mV)');
    ylabel('recovery variable, u');
    hold on;
end
