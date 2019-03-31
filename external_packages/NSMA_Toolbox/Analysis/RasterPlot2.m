function fh = RasterPlot2(S,varargin)
%  plot a SpikeRaster for a given input S-matix
%  RasterPlot2(S)
%
%  INPUT:
%    S .... S-Matrix, cell-array of ts objects
%   
%   Varargin:  Pairs of 'Parameter', Value
%       PARAMETER      VALUE
%       'Color'        any valid matlab color spec (e.g. 'r','g','k', ...)
%       'SumOfQ'  ...  0/1 logical flag
%
%
%  Version 0.1
%  Author Peter Lipa and Students of MatlabCourse 2005


%check input
if ~iscell(S) || ~strcmp(class(S{1}),'ts')
    error('RasterPlot2 Error: input S is of type %s, but must be a cell array of ts objects!', class(S));
end

% extract varargin
Extract_varargin;
if ~exist('Color','var')
    Color = 'k';
end
if ~exist('SumOfQ','var')
    SumOfQ = 1;
end


dy = 0.4;
%fh1 = figure;
clf;
if SumOfQ
   % calc Sum of Q
   Q = MakeQfromS(S,10000,'ProgressBar','none');
   Qdata = Data(Q)';
   tb = Range(Q,'sec');
   sumQ = sum(Qdata);
   ah2 = axes('position',[0.1 .7 .8 .2]);
   line(tb,sumQ);
   axis tight
end

ah1 = axes('position',[0.1 0.1 .8 .6]);
%axis(ah1);  % make ah1 the current axis

for i=1:length(S)
    t = Range(S{i});
    if size(t,2) > 1
        t = t';          % make sure t is a column vector
    end
    x = t/10000;        % plot time in sec units on x-axis
    x3 = [x,x,x];
    x3 = reshape(x3',[],1);
    
    y = i*ones(length(x),1);
    y3 = [y-dy, y+dy, y+NaN];
    y3 = reshape(y3',[],1);
    line(x3,y3,'LineStyle','-','LineWidth',0.2,'Color',Color);
    
    
end %for i

axis tight;
%title('Spike Raster');
xlabel('Time (sec)');
ylabel('Cell Number');

if nargout > 0
    fh = fh1;
end


return

