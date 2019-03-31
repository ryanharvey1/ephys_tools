% cosine_directional_tuning_preprocess

%PULL TIME STAMPS, XY COORDS AND VELOCITY FROM DATA STRUCTURE
ts=data.frames(data.frames(:,1)>data.events(1,2) & data.frames(:,1)<data.events(2,2),1);
x=data.frames(data.frames(:,1)>data.events(1,2) & data.frames(:,1)<data.events(2,2),2);
y=data.frames(data.frames(:,1)>data.events(1,2) & data.frames(:,1)<data.events(2,2),3);
vel=data.frames(data.frames(:,1)>data.events(1,2) & data.frames(:,1)<data.events(2,2),5);

%CALCULATE THE MOVMEMENT DIRECTION (BODY NOT HEAD)
[movement_angle]=XYangle(x,y);

%CREATE BINS FOR AZIMUTH
bins=[22.5:45:360-22.5,22.5];

direction={'ne_fints','n_fints', 'nw_fints', 'w_fints', 'sw_fints', 's_fints', 'se_fints','e_fints'};


for i=1:length(bins)
    if i+1>length(bins)
        break
    end
    
   idx=movement_angle>bins(i) & movement_angle<bins(i+1);
   if i==length(bins)-1
       idx=movement_angle>bins(i) | movement_angle<bins(i+1);
   end
   tempvel=vel(idx);
   tempx=x(idx);
   tempy=y(idx);
   tempts=ts(idx);
   
   idx=contiguousframes(tempvel>=7.5,12);
   tempx=tempx(idx);
   tempy=tempy(idx);
   
   [start,ends,~]=findgroups(idx);

   thetamove.dirs.(direction{i})=[tempts(start),tempts(ends)];
   thetamove.xy.(direction{i})=[tempx,tempy];
end


figure;
subplot(3,3,5)
plot(x,y,'k');hold on;box off, axis tight;axis off

place=[3 2 1 4 7 8 9 6];
for i=1:length(direction)
    subplot(3,3,place(i))
   plot(thetamove.xy.(direction{i})(:,1),thetamove.xy.(direction{i})(:,2),'.k')
end
