%PolarLabviewOutput_loop creates a polar plot from labview output .o
%
%Input:
%       - .o output file from Labview Results.table.1 cell.vi
%
%Output:
%       - figure showing polar plot of HD x FR
%
%Created by Shawn W March 2014

%identify path to data files
path = '/Users/bjclark/Desktop/Dropbox/EC_analysis/__RawData_&_Summary/__control_results_files_PaS';
Cell_list = FindFiles('*.txt', 'StartingDirectory', path);

%loops through the files and plot polar HD x FR
for i = 1:length(Cell_list);

% extract FR across bins
tFR = dlmread(Cell_list{i},'\t',10,3);
FR = transpose(tFR);
maxFR = max(FR);

% number of directional bins
dBins = 60;

% determine the size of a bin in radians
BinsAngle = (0+((2*pi/(dBins)/2)):2*pi/(dBins):(2*pi)-((2*pi/(dBins)/2)));

% Polar plot of HD x FR
subplot('Position', [0 0 1 1]);
polarplot = polar(BinsAngle([1:end 1]), FR([1:end 1]));
set(polarplot, 'linewidth',3,'color','k');
axis tight
set(0,'Showhiddenhandles','on')
extrastuff = setdiff(get(gca,'children'),polarplot);
delete(extrastuff)
horizontal=line([-100 100],[0 0]);
vertical=line([0 0],[-100 100]);
set(horizontal,'linewidth',2,'color','k');
set(vertical,'linewidth',2,'color','k');
PFR = uicontrol('Style','text','position',[585 544 90 45]);
set(PFR,'String',num2str(maxFR),'background','w','fontsize',34)
[filepath,filename] = fileparts(Cell_list{i});
saveas(polarplot,[filepath filesep filename 'polar.jpg']);
clear polarplot
clear tFR

end