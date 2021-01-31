
function HD_plots(data,session,cell,save_path,f)
% Makes a figure to check for HD firing and save to save_path destination. Uses plot functions from postprocessFigures. 
% Input: 
% - data: data structure for ephys_tools
% - session: index for spike array e.g. data.Spikes{cell,session}
% - cell: index for spike array e.g. data.Spikes{cell,session}
% - save_path: path specifying save destination. 

fig = figure; 
fig.Color = [1 1 1];
% tuning curve 
subplot(2,2,1)
postprocessFigures.plot_HD_tuning(data,session,cell)
% spike on path
subplot(2,2,3)
postprocessFigures.spikesonpath_2d(data,session,cell,'HD')
% spikes on head angle 
subplot(2,2,2)
spikesonheading(data,session,cell)
% waveform
subplot(2,2,4)
postprocessFigures.avg_waveforms(data,session,cell)

if f
    % Save the figure
    saveas(fig,[save_path,filesep,data.sessionID,'_','session_',num2str(session),'_cell_',num2str(cell),'.pdf'],'pdf')
    close all
end

end

function spikesonheading(data,session,cell)
% Grab session index 
sess_idx = data.frames(:,1) >= data.events(1,session) & data.frames(:,1) <= data.events(2,session);

% Upack angles, timestamps, and spikes
a = data.frames(sess_idx,4);
ts = data.frames(sess_idx,1);
spks = data.Spikes{cell};

% plot head angle
plot(ts,wrapTo360(a),'color',[.7,.7,.7],'LineWidth',2)
alpha(.55)
hold on
% plot spikes
scatter(spks,wrapTo360(circular_interp(ts,a,spks)),'r','filled')
alpha(.25)
legend({'Head Direction','Spike'})
ylabel('Head Angle (Degrees)')
xlabel('Time (s)')
axis tight
end