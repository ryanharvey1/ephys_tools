% MWM_IRR
clear;close all;clc
load('/Users/RyanHarvey/Downloads/NavDay5T3.mat')
KEEP('pathResults')
ID=fieldnames(pathResults);
strategylist={'Undefined','Direct','Directed','Circuitious Direct','Target Search','Target Scanning','Chaining','Foucused Search','Scanning','Scanning Surrounds','Incursion','Self-Orienting','Thigmotaxis'};
j=1;done=0;
NumofStrat=0:11;
while done==0
    ii=randi(length(ID));
    if done==1;break;end
    trial=(strsplit(ID{ii,:},'_'));
    strategytemp=pathResults.(ID{ii}).strategy;
    x=0; y=0; rad=75; th = 0:pi/179.5:2*pi;
    xunit=rad*cos(th)+x;yunit=rad*sin(th)+y;

    for i=1:length(strategytemp)
        fig=figure; fig.Color=[1 1 1]; set(fig,'Position',[1 1 960 990])
        plot(xunit,yunit,'k');hold on;axis image;box off;axis off;
        scatter(26.53, -26.53, 3000, 'k')
        plot(pathResults.(ID{ii}).segments(i).segments(:,1),pathResults.(ID{ii}).segments(i).segments(:,2));
        scatter(pathResults.(ID{ii}).segments(i).segments(1,1),pathResults.(ID{ii}).segments(i).segments(1,2), 80,'filed', 'r')
        strategy=(pathResults.(ID{ii}).strategy(i)+1);
        irr(j,1)=strategy-1;
        an=input([strategylist{strategy},'?    ']);
        if ismember(an,[1 0])==0;close all;done=1;KEEP('irr','pathResults','done','strategylist');break;end
        irr(j,2)=an;
        j=j+1;
        if length(irr)>length(ID)*.15;disp('You have labeled 15% of data');end
        close all
    end
end
clc
disp([num2str(length(irr(:,2))),' Labeled'])
disp([num2str(sum(irr(:,2))),' Yes'])
disp([num2str(length(irr(irr(:,2)==0,1))),' No'])
disp(['Score: ',num2str(sum(irr(:,2))/length(irr(:,2)))])
disp(['Most Incorrectly Classified: ',strategylist{mode(irr(irr(:,2)==0,1))}])

function KEEP(varargin)
if isempty(varargin);return;end;wh = evalin('caller','who');
if isempty(wh);error('  There is nothing to keep!');end
variable = [];for i = 1:length(wh);variable = [variable,':',wh{i}];end;variable = [variable,':'];flag = 0;
for i = 1:length(varargin)
    I = findstr(variable,[':',varargin{i},':']);
    if isempty(I)
        flag = 1;
    elseif I == 1
        variable = variable(1+length(varargin{i})+1:length(variable));
    elseif I+length(varargin{i})+1 == length(variable)
        variable = variable(1:I);
    else
        variable = [variable(1:I),variable(I+length(varargin{i})+2:length(variable))];
    end
end
if flag == 1;return;end;I = findstr(variable,':');
if length(I) ~= 1
    for i = 1:length(I)-1
        if i ~= length(I)-1;del(i) = {[variable(I(i)+1:I(i+1)-1),' ']};else;del(i) = {variable(I(i)+1:length(variable)-1)};end
    end
    evalin('caller',['clear ',del{:}])
end
end