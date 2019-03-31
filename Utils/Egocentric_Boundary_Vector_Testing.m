% Egocentric_Boundary_Vector_Testing
com=which('Egocentric_Boundary_Vector_Testing');
com=strsplit(com,filesep);
basedir=[com{1},filesep,'Users',filesep,com{3},filesep,'GoogleDrive',filesep,'MatlabDir'];
addpath([basedir,filesep,'polarPcolor'])

%% Manual setings 
% set data path
data=load('D:\Projects\HPCatn\ProcessedData\HPCatn05_S20181115110624');

session=2;

% set n directional bins (6 degree bins)
thetabins=61;
% set n distance bins (3cm bins)
distancebins = (round(data.maze_size_cm(session)/3)+1);

%%
theta=0:.01:2*pi;
c=hsv(length(theta));



ts=data.frames(data.frames(:,1)>data.events(1,session) & data.frames(:,1)<data.events(2,session),1);
x=data.frames(data.frames(:,1)>data.events(1,session) & data.frames(:,1)<data.events(2,session),2);
y=data.frames(data.frames(:,1)>data.events(1,session) & data.frames(:,1)<data.events(2,session),3);
vel=data.frames(data.frames(:,1)>data.events(1,session) & data.frames(:,1)<data.events(2,session),5);

ts(isnan(x))=[];
x(isnan(x))=[];
y(isnan(y))=[];
vel(isnan(vel))=[];

in=contiguousframes(vel<5,1);


K = convhull(x,y);

xbound=x(K);
ybound=y(K);




xbound=interp1(linspace(1,1000,length(xbound)),xbound,1:1000);
ybound=interp1(linspace(1,1000,length(ybound)),ybound,1:1000);
% calculate closest distance to boundary for all points
distances=zeros(1,length(x));
for i=1:length(x)
    distances(i)=min(sqrt(sum(bsxfun(@minus, [xbound;ybound]', [x(i),y(i)]).^2,2)));
end

% find movement heading
[ angle ] = XYangle(x,y);
angle=[angle;angle(end)];
% no movement will give you 0 or 180 degrees
angle(angle==0 | angle==180)=NaN;
% smooth over no movement 
angle=wrapTo360(fixNLXangle(angle,round(0.1667*data.samplerate)));

% filter out low velocity 
ts(in)=[];
x(in)=[];
y(in)=[];
vel(in)=[];
angle(in)=[];
distances(in)=[];


% bin occ by theta and distance
occ = histcounts2(angle',distances',linspace(0,360,thetabins),linspace(min(distances),max(distances),distancebins));
filtWidth = [7,7]; filtSigma = 5;
imageFilter=fspecial('gaussian',filtWidth,filtSigma);
occ=conv2([occ;occ;occ],imageFilter,'same')/data.samplerate;
occ=occ(thetabins:(thetabins-1)*2,:);


ts=data.frames(data.frames(:,1)>data.events(1,session) & data.frames(:,1)<data.events(2,session),1);
x=data.frames(data.frames(:,1)>data.events(1,session) & data.frames(:,1)<data.events(2,session),2);
y=data.frames(data.frames(:,1)>data.events(1,session) & data.frames(:,1)<data.events(2,session),3);
vel=data.frames(data.frames(:,1)>data.events(1,session) & data.frames(:,1)<data.events(2,session),5);
vel(isnan(x))=[];
ts(isnan(x))=[];
x(isnan(x))=[];
y(isnan(y))=[];

in=contiguousframes(vel<5,1);

[ angle ] = XYangle(x,y);
angle=[angle;angle(end)];
angle(angle==0 | angle==180)=NaN;
angle=wrapTo360(fixNLXangle(angle,round(0.1667*data.samplerate)));

distances=zeros(1,length(x));
for i=1:length(x)
    distances(i)=min(sqrt(sum(bsxfun(@minus, [xbound;ybound]', [x(i),y(i)]).^2,2)));
end

[p,n]=numSubplots(length(data.Spikes));
for i=1:length(data.Spikes)
% cell to look at
spkts=data.Spikes{i};  
spkts=spkts(spkts>data.events(1,session) & spkts<data.events(2,session));


% bin spikes by theta and distance
anglespk=interp1(ts,angle,spkts);
distancesspk=interp1(ts,distances,spkts);
inspk=interp1(ts,double(in),spkts);

anglespk(inspk>0)=[];
distancesspk(inspk>0)=[];
spkts(inspk>0)=[];


spkmap = histcounts2(anglespk,distancesspk,linspace(0,360,thetabins),...
    linspace(min(distances),max(distances),distancebins));
imageFilter=fspecial('gaussian',filtWidth,filtSigma);
spkmap=conv2([spkmap;spkmap;spkmap],imageFilter,'same');
spkmap=spkmap(thetabins:(thetabins-1)*2,:);

% divide to get firing rate
map=(spkmap./occ);
map(isnan(map))=0;
map(isinf(map))=0;
% map=gpuArray(map);

% %% data plot
% figure
% R = linspace(0,max(distances),size(map,2)); % (distance in km)
% theta = round(linspace(0,360,size(map,1))); % in degrees
% Z=map';
% 
% rMin = min(R);
% rMax = max(R);
% thetaMin=min(theta);
% thetaMax =max(theta);
% % Definition of the mesh
% Rrange = rMax - rMin; % get the range for the radius
% rNorm = R/Rrange; %normalized radius [0,1]
% % get hold state
% cax = newplot;
% % transform data in polar coordinates to Cartesian coordinates.
% YY = (rNorm)'*cosd(theta);
% XX = (rNorm)'*sind(theta);
% % plot data on top of grid
% h = pcolor(XX,YY,Z,'parent',cax);
%  shading flat
% colormap jet
% axis image
% shading interp
% 
% figure
% R = linspace(0,max(distances),size(map,2)); 
% Az = linspace(0,360,size(map,1)); % in degrees
% [rg, thg] = meshgrid(R,deg2rad(round(Az)));
% [x,y] = pol2cart(thg,rg);
% pcolor(y,x,map);
% colormap jet;
% shading interp;
% axis tight;
% axis image

map=imresize(map,[size(map)*20]);


% figure;
% polarscatter(deg2rad(anglespk),distancesspk,10,deg2rad(anglespk));
% colormap(c)


% plot using polarPcolor  
figure(1)
subplot(p(1),p(2),i);
R = linspace(0,max(distances),size(map,2)); 
% R = logspace(log10(0.1),log10(max(distances)),size(map,2)); 
Az = linspace(0,360,size(map,1)); % in degrees
% [~,~]= polarPcolor(R,round(Az),map,'Rscale','log','Nspokes',5,'Ncircles',2);
[~,~]= polarPcolor(R,round(Az),map,'Nspokes',5,'Ncircles',2);
colorbar off
colormap jet
% shading interp
axis tight


figure(2)
subplot(p(1),p(2),i);
plot(x,y,'.k');
hold on;
scatter(interp1(ts,x,spkts),interp1(ts,y,spkts),10,...
   anglespk,'filled');
ax=gca;
colormap(ax,c)
box off
axis off
axis image
end


% N = 360;
% R = linspace(0,1000,N)./1000; % (distance in km)
% Az = linspace(0,360,N); % in degrees
% [~,~,windSpeed] = peaks(N); % radial wind speed
% figure
% [~,c]= polarPcolor(R,Az,windSpeed);
% ylabel(c,' radial wind speed (m/s)');
% set(gcf,'color','w')
% 
% [h,c]=polarPcolor(distances,angle,windSpeed)



