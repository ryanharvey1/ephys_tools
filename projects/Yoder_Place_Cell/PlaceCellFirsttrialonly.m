% PlaceCellFirsttrialonly
clear;clc;close all
A = importdata('/Users/ryanharvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/Place cell First trial only .xlsx') ;
%%
mice=fieldnames(A.data);

controlmice=mice(contains(mice,'Pass'));
control=zeros(1,5);
for i=1:length(controlmice)
    temp=A.data.(controlmice{i});
    control=[control;temp(:,[5:8,10])];
    
end
control(1,:)=[];
control(control>10000)=control(control>10000)/1000000;
control(control<0)=abs(control(control<0));


tiltedmice=mice(contains(mice,'Fail'));
tilted=zeros(1,5);
for i=1:length(tiltedmice)
    temp=A.data.(tiltedmice{i});
    tilted=[tilted;temp(:,[5:8,10])];
    
end
tilted(1,:)=[];
tilted(tilted>10000)=tilted(tilted>10000)/1000000;
tilted(tilted(:,2)>1,2)=tilted(tilted(:,2)>1,2)/100;
tilted(tilted<0)=abs(tilted(tilted<0));

% Need to locate AP durations for missing units
% tilted(tilted(:,5)==0,4)
% tointerp=tilted(:,4:5);
% tointerp(logical(sum(isnan(tointerp),2)),:)=[];
% tointerp(logical(sum(tointerp==0,2)),:)=[];
% interp1(tointerp(:,1),tointerp(:,2),tilted(tilted(:,5)==0,4),'nearest')
% figure;scatter(tointerp(:,1),tointerp(:,2))

%% Kmeans to pull out principle cells 
% idx = kmeans(tilted,3);
% 
% figure
% scatter(tilted(idx==1,4),tilted(idx==1,5));hold on
% scatter(tilted(idx==2,4),tilted(idx==2,5))
% scatter(tilted(idx==3,4),tilted(idx==3,5))
% 
% figure
% silhouette(tilted,idx)
% 
% %
% idx = kmeans(control,3);
% 
% figure
% scatter(control(idx==1,4),control(idx==1,5));hold on
% scatter(control(idx==2,4),control(idx==2,5))
% scatter(control(idx==3,4),control(idx==3,5))
% 
% figure
% silhouette(control,idx)


%% FIND GOOD CRITERIA
% control2=control(control(:,4)<10 & control(:,5)>185,:);
% tilted2=tilted(tilted(:,4)<10 & tilted(:,5)>185,:);
% 
% ncells=sum([tilted2(:,2)>.5 & tilted2(:,3)>.6;control2(:,2)>.5 & control2(:,3)>.6]);
% disp(['Number of Cells: ',num2str(ncells)])

%% SEPERATE PRINCIPLE CELLS
% count all cells
nallcells_control=size(control,1);
nallcells_tilted=size(tilted,1);

% find cells that fire <=10 hz on average (cells that fire quickly and continuously tend to be interneurons)
control=control(control(:,4)<=10,:);
tilted=tilted(tilted(:,4)<=10,:);

% find cells with an AP duration of >= 200 microseconds (cells with AP durations less than 200 us tend to be fast spiking interneurons)
% The use of the or opperator is used to be liberal with NaNs and 0s 
control=control(control(:,5)>=185 | isnan(control(:,5)) | control(:,5)==0,:);
tilted=tilted(tilted(:,5)>=185 | isnan(tilted(:,5)) | tilted(:,5)==0,:);

% count principle cells
nprinciplecells_control=size(control,1);
nprinciplecells_tilted=size(tilted,1);

[h,p, chi2stat,df] = prop_test([nprinciplecells_control nprinciplecells_tilted] , [nallcells_control nallcells_tilted], 0)
disp(['% of control principle cell: ',num2str(nprinciplecells_control/nallcells_control)])
disp(['% of tilted principle cell: ',num2str(nprinciplecells_tilted/nallcells_tilted)])

[h,p, chi2stat,df] = prop_test([54 100] , [nprinciplecells_control nprinciplecells_tilted], 0)
disp(['% of control place cells: ',num2str(54/nprinciplecells_control)])
disp(['% of tilted place cells: ',num2str(100/nprinciplecells_tilted)])


%% PLOT PRINCIPLE CELLS
close all
load('/Users/ryanharvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/PlaceCellFeatures4PlaceCellFirstTrialOnly.mat')

controlexamples=[3,5,6,11,12,17,31];
tiltedexamples=[19,38,43,68,69,75,77];
varnames={'Peak Rate (hz)','Coherence','Information Content (bits/spk)','Average Rate (hz)'};
for i=1:4   
    figure(i)
    [f1,x1] = ecdf(control(:,i));
    [f2,x2] = ecdf(tilted(:,i));
    
    p1=plot(x1,f1);
    set(p1,'LineWidth',4,'Color','k')
    hold on
    examplecontrol=controlPlaceCells(controlexamples,i);
    
    for ii=1:7
        plot([examplecontrol(ii);examplecontrol(ii)],[0;1],'k');hold on
    end
%     examplecontrol(ii)==f1
    
    
    p2=plot(x2,f2);
    set(p2,'LineWidth',4,'Color','r')
    exampletilt=tiltedPlaceCells(tiltedexamples,i);
    
    for ii=1:7
        plot([exampletilt(ii);exampletilt(ii)],[0;1],'r');hold on
    end        
    
    ylabel('Cumulative Frequency')
    xlabel(varnames{i})
    
    set(gca,'FontSize',20,'FontWeight','bold','LineWidth',2,'box','off')
    

    
end


%% scatter plots of distribution
load('/Users/ryanharvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/PlaceCellFeatures4PlaceCellFirstTrialOnly.mat')
close all

% plot all cells
scatter_fig=figure; scatter_fig.Color=[1 1 1];
plot(control(:,3),control(:,2),'.k');hold on
plot(tilted(:,3),tilted(:,2),'.r');
xlabel('Information Content (bits/spk)')
ylabel('Coherence')
set(gca,'Box','off','FontWeight','bold','FontSize',15,'LineWidth',3)

% plot cut off lines
% l1=plot([.6;.6],[0;1],'k');
% l2=plot([0;max([tilted(:,3);control(:,3)])],[.5;.5],'k');
% set(l1,'LineWidth',3)
% set(l2,'LineWidth',3)

% plot place cells
scatter(controlPlaceCells(:,3),controlPlaceCells(:,2),40,'k','filled');hold on
scatter(tiltedPlaceCells(:,3),tiltedPlaceCells(:,2),40,'r','filled');


% arrow to example cells
controlexamples=[3,5,6,11,12,17,31];
tiltedexamples=[19,38,43,68,69,75,77];

% dim = [.2 .5 .3 .3];
% 
% a=annotation('textbox',dim,'String','1','FitBoxToText','on');
% 
% 
% for i=1:length(controlexamples)
%     l1=plot([controlPlaceCells(controlexamples(i),3);controlPlaceCells(controlexamples(i),3)],[controlPlaceCells(controlexamples(i),2),1],'k')
%     set(l1,'LineWidth',3)
%     
%     dim=[controlPlaceCells(controlexamples(i),3)/10, .5, .3, .3]
%     a=annotation('textbox',dim,'String',num2str(i),'FitBoxToText','on');
% 
% end
% 
% for i=1:length(tiltedexamples)
%     l1=plot([tiltedPlaceCells(tiltedexamples(i),3);tiltedPlaceCells(tiltedexamples(i),3)],[tiltedPlaceCells(tiltedexamples(i),2),1],'r')
%     set(l1,'LineWidth',3)
% end

% temp=[rescale([controlPlaceCells(:,2);tiltedPlaceCells(:,2)],0,1),rescale([controlPlaceCells(:,3);tiltedPlaceCells(:,3)],0,1)];

% control4arrow=temp(1:size(controlPlaceCells),:);
% tilted4arrow=temp(1:size(tiltedPlaceCells),:);
% 
% x=[controlPlaceCells(3,3)]
% y=[controlPlaceCells(3,2)]
% 
% x=[.8,.7];
% y=[.8,.7];
% % 
% annotation('textarrow',x,y,'String','Cell 2');




%%
% scatter(tilted(:,4),tilted(:,5))

% load('/Users/ryanharvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/Place Cell First trial only.mat')


%% Peak Rate
figure;
[f,x]=ecdf(control.peakrate);

plot(x,f);hold on

[f,x]=ecdf(tilted.peakrate);

plot(x,f);
xlabel('Peak Rate')
ylabel('Cumulative Frequency')

%% coherence
figure;
[f,x]=ecdf(control.coherence);

plot(x,f);hold on

[f,x]=ecdf(tilted.coherence);

plot(x,f);
xlabel('Coherence')
ylabel('Cumulative Frequency')
%% infocontent
figure;
[f,x]=ecdf(control.infocontent);

plot(x,f);hold on

[f,x]=ecdf(tilted.infocontent);

plot(x,f);
xlabel('Information Content')
ylabel('Cumulative Frequency')
%% avgrate
figure;
[f,x]=ecdf(control.avgrate);

plot(x,f);hold on

[f,x]=ecdf(tilted.avgrate);

plot(x,f);
xlabel('Average Rate')
ylabel('Cumulative Frequency')
%% waveformdur
figure;
[f,x]=ecdf(control.waveformdur);

plot(x,f);hold on

[f,x]=ecdf(tilted.waveformdur);

plot(x,f);
xlabel('Wave Form Duration')
ylabel('Cumulative Frequency')
%%



% createcdf(control)
% createcdf(tilted)


% function createcdf(data)
% measures=fieldnames(data);
% for i=1:length(measures)
%     data.(measures{i});
% end
% 
% 
% end