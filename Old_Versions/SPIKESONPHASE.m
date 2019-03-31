% SPIKESONPHASE
clc, clear, close all
tableexample = readtable('F:\Users\reharvey\Place_Cell_Data\PAE_Rat\_AllSpikeDataRH11.xlsx');
tableexample(:,1)=[];
% ###################
path='F:\Users\reharvey\Place_Cell_Data\PAE_Rat';
RH11 = xlsread('F:\Users\reharvey\Place_Cell_Data\PAE_Rat\_AllSpikeDataRH11.xlsx');
RH16 = xlsread('F:\Users\reharvey\Place_Cell_Data\PAE_Rat\_AllSpikeDataRH16.xlsx');
RH13 = xlsread('F:\Users\reharvey\Place_Cell_Data\PAE_Rat\_AllSpikeDataRH13.xlsx');
RH14 = xlsread('F:\Users\reharvey\Place_Cell_Data\PAE_Rat\_AllSpikeDataRH14.xlsx');

control=[RH13;RH14];
PAE=[RH11;RH16];

% filter spikes
control=control(control(:,15)>=50,:);
PAE=PAE(PAE(:,15)>=50,:);

% Filter by information content
control=control(control(:,1)>=.30,:);
PAE=PAE(PAE(:,1)>=.30,:);
%

ii=1;
variables={'DELTA','THETA','ALPHA','BETA','GAMMA','HI-GAMMA'};
for iii=1:2
    if iii==1;group=control; elseif iii==2;group=PAE; end;
    ii=1;
    for i=40:45
%         rayleighsig=28:33;
%         group=group(group(:,rayleighsig(ii))>=0.05,:);
        if iii==1; figure(i); elseif iii==2; figure(i+6); end;
        % plot wave
        steps=(2*pi)/(length(group)-1);
        t = 0:steps:(2*pi);
        y = sin(t);
        hh=plot(t,y,'k'); hold on
        hh.LineWidth=5;
        set(gca, 'XTick', []);
        set(gca, 'YTick', []);
        box off
        title(variables(ii),'FontSize',12); ii=ii+1;
        % scatter spikes
        spikes=sin(deg2rad(group(:,i))');
        h=scatter(deg2rad(group(:,i))',spikes,'Marker','x');
        h.MarkerEdgeColor=[1 0 0];
        h.MarkerFaceColor=[1 0 0];
        h.LineWidth = 8;
        axis([0 2*pi -1.2 1.2])
    end
end

PAEperTheta=(sum(PAE(:,29)<=0.05)/length(PAE))*100;
contrlperTheta=(sum(control(:,29)<=0.05)/length(control))*100;

PAEperGamma=(sum(PAE(:,32)<=0.05)/length(PAE))*100;
contrlperGamma=(sum(control(:,32)<=0.05)/length(control))*100;

PAEperHIGamma=(sum(PAE(:,33)<=0.05)/length(PAE))*100;
contrlperHIGamma=(sum(control(:,33)<=0.05)/length(control))*100;
