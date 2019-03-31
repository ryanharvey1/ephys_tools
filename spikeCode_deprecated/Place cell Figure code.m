%% Matlab Code For Place cell Figures

%Code for Avg Firing Rate

[d,e]=hist(controlfiltii(:,1),50);
 sumd=sum(d); dnorm=(d/sumd)*100;
 figure (1)
 plot(e,dnorm,'k')

[f,g]=hist(tiltfiltii(:,1),50);
sumf=sum(f); fnorm=(f/sumf)*100;
figure (2)
plot(g,fnorm,'k')

figure (1),hold on
 plot (g,fnorm,'r')

%% Code for For Spatial Coherence graph

[d,e]=hist(controlfilti(:,2),50);
 sumd=sum(d); dnorm=(d/sumd)*100;
 figure (5)
 plot(e,dnorm,'k')
hold on
xlim([0 1])

[f,g]=hist(tiltfilti(:,2),50);
sumf=sum(f); fnorm=(f/sumf)*100;
figure (6)
plot(g,fnorm,'k')
hold on
 xlim([0 1])

figure (5),hold on
 plot (g,fnorm,'r')


%% Code for Spatial Information Content

[d,e]=hist(controlfiltii(:,3),50);
 sumd=sum(d); dnorm=(d/sumd)*100;
 figure (5)
 plot(e,dnorm,'k')

[f,g]=hist(tiltfiltii(:,3),50);
sumf=sum(f); fnorm=(f/sumf)*100;
figure (6)
plot(g,fnorm,'k')

figure (5),hold on
 plot (g,fnorm,'r')


%% Code for Avg Rate

[d,e]=hist(controlfilti(:,4),50);
 sumd=sum(d); dnorm=(d/sumd)*100;
 figure (5)
 plot(e,dnorm,'k')

[f,g]=hist(tiltfilti(:,4),50);
sumf=sum(f); fnorm=(f/sumf)*100;
figure (6)
plot(g,fnorm,'k')

figure (5),hold on
 plot (g,fnorm,'r')

%% Code for Scatter Plot

figure (3)
scatter(tiltfilti(:,4),tiltfilti(:,3))
figure (3), hold on
lsline

%% Code for greater than / less than filters

x=find(controlfilti(:,4)>.4); controlfiltii=controlfilti(x,:); 
y=find(tiltfilti(:,4)>.4); tiltfiltii=tiltfilti(y,:);