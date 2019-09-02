
%standard error of linear and angular velocity
param_idx=params.subID;

for i=1:length(param_idx)
    runIdx=params.runIdx{i};
    temp_angVel=params.angVel{i};
    temp_linVel=params.pathIV{i};
    
    if contains(param_idx(i),'Tg')
        motion(i,1)=1;
    else
        motion(i,1)=0;
    end
    
    if mod(i,2)==0
        motion(i,2)=2;
    else
        motion(i,2)=1;
    end
    
    motion(i,3)=nanstd(abs(temp_angVel));
    motion(i,4)=nanstd(abs(temp_linVel));
    motion(i,5)=nanstd(abs(temp_angVel(runIdx)));
    motion(i,6)=nanstd(abs(temp_linVel(runIdx)));
end

figure; 
subplot(1,2,1)
plotspread_wrapper(motion(motion(:,1)==0 & motion(:,2)==1,3),motion(motion(:,1)==1 & motion(:,2)==1,3),{'WT','Tg'})
title('Standard Deviation Angular Velocity - Day 1')
ylabel('deg/s')
ylim([0 max(motion(:,5))])
subplot(1,2,2)
plotspread_wrapper(motion(motion(:,1)==0 & motion(:,2)==2,3),motion(motion(:,1)==1 & motion(:,2)==2,3),{'WT','Tg'})
title('Standard Deviation Angular Velocity- Day 2')
ylabel('deg/s')
ylim([0 max(motion(:,5))])

figure; 
subplot(1,2,1)
plotspread_wrapper(motion(motion(:,1)==0 & motion(:,2)==2,4),motion(motion(:,1)==1& motion(:,2)==2,4),{'WT','Tg'})
title('Standard Deviation linear Velocity - Day 1')
ylabel('deg/s')
ylim([0 max(motion(:,6))])
subplot(1,2,2)
plotspread_wrapper(motion(motion(:,1)==0 & motion(:,2)==1,4),motion(motion(:,1)==1 & motion(:,2)==1,4),{'WT','Tg'})
title('Standard Deviation linear Velocity- Day 2')
ylabel('deg/s')
ylim([0 max(motion(:,6))])



stat_plot(motion(motion(:,1)==0 & motion(:,2)==2,3),motion(motion(:,1)==1 & motion(:,2)==2,3),{'WT','Tg'},'Angular Velocity','plots',2)

stat_plot(motion(motion(:,1)==0 & motion(:,2)==2,4),motion(motion(:,1)==1 & motion(:,2)==2,4),{'WT','Tg'},'Linear Velocity','plots',2)


